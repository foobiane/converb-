#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "m_pd.h"
#include "fftw3.h"

#define HOP_FAC 0.5

/*
 * A note on block size
 *
 * Block size can change during the execution of this external. This is treated
 * as an error, and is handled accordingly. A real-time solution could be
 * devised, though it would effectively involve reloading the impulse response
 * during DSP. This can cause massive DSP slowdowns, and there's minimal
 * performance benefit compared to stopping DSP and reloading from scratch.
 * 
 */

typedef struct _converb {
    t_object obj;

    float* input_buffer;
    fftwf_complex* ir_dft;
    float* circular_output_buffer;

    int output_buffer_size;
    int output_read;
    int output_write;

    int num_samples_read;

    int data_size; // num samples for all DFTs, sample arrays, and buffers
    int block_size; 
    int hop_size;

    bool loaded;

    float fsig;
} t_converb;

static t_class* converb_class;

// Pd Prototypes
void converb_tilde_setup(void);
void* converb_new(t_symbol* s, short argc, t_atom* argv);
void converb_free(t_converb* x);
void converb_load(t_converb* x, t_symbol* s);
void converb_dsp(t_converb* x, t_signal** sp);
t_int* converb_perform(t_int* w);

// Helper Functions
void converb_write_to_output_buffer(t_converb* x, float* window);
void converb_read_from_output_buffer(t_converb* x, float* buf, int n);
float* convolve_with(float* data, fftwf_complex* ir_dft, int n);
float* generate_hann_window(int n);
void normalize(float* x, int size);
void window(float* data, int n);

static float* hann_window = NULL;

///////////////////////////////////////////////////////////////////////////////
// PD FUNCTIONALITY
///////////////////////////////////////////////////////////////////////////////

void converb_tilde_setup(void) {
    converb_class = class_new(
        gensym("converb~"),
        (t_newmethod) converb_new,
        (t_method) converb_free,
        sizeof(t_converb),
        CLASS_DEFAULT,
        A_GIMME,
        0
    );

    CLASS_MAINSIGNALIN(converb_class, t_converb, fsig);

    class_addmethod(converb_class, (t_method) converb_dsp, gensym("dsp"), A_CANT, 0);	
    class_addmethod(converb_class, (t_method) converb_load, gensym("load"), A_SYMBOL, 0);
}

void* converb_new(t_symbol* s, short argc, t_atom* argv) {
    t_converb* x = (t_converb*) pd_new(converb_class);
    outlet_new(&x->obj, gensym("signal"));

    // Setting all fields to their default values
    x->input_buffer = NULL;
    x->circular_output_buffer = NULL;
    x->ir_dft = NULL;
    
    x->output_buffer_size = 0;
    x->output_read = 0;
    x->output_write = 0;

    x->num_samples_read = 0;

    x->data_size = 0;
    x->block_size = sys_getblksize(); // default block size is the system's
    x->hop_size = 0;

    x->loaded = false;

    post("converb~: Initialized with block_size = %d", x->block_size); // debug

    return x;
}

void converb_free(t_converb* x) {
    if (x->input_buffer) free(x->input_buffer);
    if (x->circular_output_buffer) free(x->circular_output_buffer);
    if (x->ir_dft) fftwf_free(x->ir_dft);
}

void converb_load(t_converb* x, t_symbol* s) {
    post("converb~: Loading IR...");
    t_garray* arr = (t_garray*) pd_findbyclass(s, garray_class);

    if (!arr) {
        pd_error(x, "converb~: Array not found.");
        return;
    }

    int num_samples;
    t_word* sample_ptr;

    garray_getfloatwords(arr, &num_samples, &sample_ptr);
    post("converb~: num_samples = %d, sample_ptr = %p", num_samples, sample_ptr); // debug

    // Padding zeros to make the IR a multiple of the block size
    int ir_size = ceil((float) num_samples / (float) x->block_size) * x->block_size;
    float* ir = calloc(ir_size, sizeof(float));

    for (int i = 0; i < ir_size; i++)
        if (i < num_samples)
            ir[i] = (float) sample_ptr[i].w_float;
        else
            ir[i] = 0.0;

    // Initializing
    x->data_size = ir_size;
    x->hop_size = ceil((float) (HOP_FAC * ir_size) / (float) x->block_size) * x->block_size;
    x->output_buffer_size = x->data_size + x->hop_size;

    // Windowing stuff
    if (hann_window) free(hann_window);
    hann_window = generate_hann_window(ir_size);
    window(ir, ir_size);

    // Initializing buffers
    if (x->input_buffer) free(x->input_buffer);
    x->input_buffer = calloc(x->data_size, sizeof(float));

    if (x->circular_output_buffer) free(x->circular_output_buffer);
    x->circular_output_buffer = calloc(x->output_buffer_size, sizeof(float));

    if (x->ir_dft) fftwf_free(x->ir_dft);
    x->ir_dft = (fftwf_complex*) fftwf_alloc_complex(x->data_size);

    // Taking DFT of the IR and storing it in x
    fftwf_plan p = fftwf_plan_dft_r2c_1d(x->data_size, ir, x->ir_dft, FFTW_ESTIMATE);
    fftwf_execute(p);

    // Cleanup
    fftwf_destroy_plan(p);
    free(ir);

    post("converb~: Impulse response DFT loaded."); // debug
    x->loaded = true;
}

void converb_dsp(t_converb* x, t_signal** sp) {
    if (!x->ir_dft) {
        pd_error(x, "converb~: Impulse response not loaded. DSP skipped.");
        x->loaded = false; // sets flag to avoid repeating prints
    }
    else if (!x->loaded)
        return;
    else
        dsp_add(converb_perform, 4, x, sp[0]->s_vec /* Left inlet */, sp[1]->s_vec /* Outlet */, sp[1]->s_n /* Signal vector size */); 
}

t_int* converb_perform(t_int* w) {
    t_converb* x = (t_converb*) w[1];
    float* in = (float*) w[2];
    float* out = (float*) w[3];
    int n = (int) w[4];

    // Handling changing block size during DSP; note that this behavior is considered to be an error.
    if (x->block_size != n) {
        pd_error(x, "converb~: Block size changed during DSP. Please end DSP and reload the impulse response.");
        x->block_size = n; // setting our block size to the new one
        x->loaded = false; // mark IR as unloaded
        
        memset(out, 0.0, n * sizeof(float));
        return w + 5; // do nothing
    }
    else if (!x->loaded) {
        memset(out, 0.0, n * sizeof(float));
        return w + 5; // do nothing
    }

    // READ STEP: We always read from the input
    memcpy(&x->input_buffer[x->num_samples_read], in, n * sizeof(float));
    x->num_samples_read += n;

    // CONVOLUTION: If our input buffer is full, convolve it with the IR and reset the counter
    if (x->num_samples_read == x->data_size) {
        float* convolution = convolve_with(x->input_buffer, x->ir_dft, x->data_size);
        converb_write_to_output_buffer(x, convolution);

        int amount_overlap = x->data_size - x->hop_size;
        memcpy(x->input_buffer, &x->input_buffer[x->hop_size], amount_overlap * sizeof(float));

        x->num_samples_read = amount_overlap;
    }
    
    // WRITE STEP: We always write to the output
    converb_read_from_output_buffer(x, out, n);

    return w + 5;
}

void converb_write_to_output_buffer(t_converb* x, float* window) {
    for (int i = 0; i < x->data_size; i++) {
        if (i < x->data_size - x->hop_size)
            x->circular_output_buffer[(x->output_write + i) % x->output_buffer_size] += window[i];
        else
            x->circular_output_buffer[(x->output_write + i) % x->output_buffer_size] = window[i];
    }

    x->output_write = (x->output_write + x->hop_size) % x->output_buffer_size;
}

void converb_read_from_output_buffer(t_converb* x, float* buf, int n) {
    for (int i = 0; i < n; i++)
        buf[i] = x->circular_output_buffer[(x->output_read + i) % x->output_buffer_size];

    x->output_read = (x->output_read + n) % x->output_buffer_size;
}

float* convolve_with(float* data, fftwf_complex* ir_dft, int n) {
    float* copy = calloc(n, sizeof(float));
    memcpy(copy, data, n * sizeof(float));
    data = copy;

    window(data, n);

    // Taking DFT of input; if no input DFT exists yet, make space
    fftwf_complex* data_dft = fftwf_alloc_complex(n);
    
    // TODO: Potential optimization of saving this plan in the converb~ object
    fftwf_plan p1 = fftwf_plan_dft_r2c_1d(n, data, data_dft, FFTW_ESTIMATE);
    fftwf_execute(p1);
    fftwf_destroy_plan(p1);

    free(copy);

    // Point-wise multiplication of IR and input DFTs
    // Note that the input DFT stores the results of this multiplication
    // A potential optimization could be the use of SIMD intrinsics to shave off some factors of n, but I have no clue how to do that!
    for (int i = 1; i < n / 2 + 1; i++)
        data_dft[i] *= ir_dft[i];
        
    // Taking inverse DFT of the convolution; if no convolution has been done yet, make space
    float* result = calloc(n, sizeof(float));

    // Inverse DFT
    // TODO: Potential optimization of saving this plan in the converb~ object
    fftwf_plan p2 = fftwf_plan_dft_c2r_1d(n, data_dft, result, FFTW_ESTIMATE);
    fftwf_execute(p2);
    fftwf_destroy_plan(p2);

    window(result, n);

    return result;

    // normalize(x->convolution, x->data_size);
}

void window(float* data, int n) {
    for (int i = 0; i < n; i++)
        data[i] *= hann_window[i];
}

// Generates a Hann window for windowing our samples.
float* generate_hann_window(int n) {
    float* result = calloc(n, sizeof(float));

    for (int i = 0; i < n; i++)
        result[i] = sqrt(pow(sin(M_PI * i / n), 2));

    return result;
}

void normalize(float* x, int size) {
    float max = 0;

    for (int i = 0; i < size; i++) {
        if (fabs(x[i]) > max)
            max = fabs(x[i]);
    }

    if (max != 0.0) {
        for (int i = 0; i < size; i++)
            x[i] /= (max * 2);
    }
}