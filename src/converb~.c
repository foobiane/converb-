#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "m_pd.h"
#include "fftw3.h"

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
    float* convolution;
    fftwf_complex* input_dft;
    fftwf_complex* ir_dft;

    int num_samples_read;
    int num_samples_written;

    int data_size; // num samples for all DFTs, sample arrays, and buffers
    int block_size; 

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
void converb_convolve(t_converb* x);
void normalize(float* x, int size);
void complex_multiply(float c1[2], float c2[2]);

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
    x->convolution = NULL;
    x->input_dft = NULL;
    x->ir_dft = NULL;
    
    x->num_samples_read = 0;
    x->num_samples_written = 0;

    x->data_size = 0;
    x->block_size = sys_getblksize(); // default block size is the system's

    x->loaded = false;

    post("converb~: Initialized with block_size = %d", x->block_size); // debug

    return x;
}

void converb_free(t_converb* x) {
    if (x->input_buffer) free(x->input_buffer);
    if (x->convolution) free(x->convolution);
    if (x->input_dft) fftwf_free(x->input_dft);
    if (x->ir_dft) fftwf_free(x->input_dft);

    // More?
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

    if (x->input_buffer) free(x->input_buffer);
    x->input_buffer = calloc(x->data_size, sizeof(float));

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
        converb_convolve(x);
        x->num_samples_read = 0;
    }
    
    // WRITE STEP: If we have convolved data to write, copy it to the output; else, set the out to zeros.
    if (!x->convolution)
        memset(out, 0.0, n * sizeof(float));
    else {
        memcpy(out, &x->convolution[x->num_samples_written], n * sizeof(float));
        x->num_samples_written = (x->num_samples_written + n) % x->data_size;
    }

    return w + 5;
}

void converb_convolve(t_converb* x) {
    // Taking DFT of input; if no input DFT exists yet, make space
    if (!x->input_dft)
        x->input_dft = fftwf_alloc_complex(x->data_size);
    
    // TODO: Potential optimization of saving this plan in the converb~ object
    fftwf_plan p1 = fftwf_plan_dft_r2c_1d(x->data_size, x->input_buffer, x->input_dft, FFTW_ESTIMATE);
    fftwf_execute(p1);
    fftwf_destroy_plan(p1);

    // Point-wise multiplication of IR and input DFTs
    // Note that the input DFT stores the results of this multiplication
    // A potential optimization could be the use of SIMD intrinsics to shave off some factors of n, but I have no clue how to do that!
    for (int i = 0; i < x->data_size / 2 + 1; i++)
        complex_multiply(x->input_dft[i], x->ir_dft[i]);

    // Taking inverse DFT of the convolution; if no convolution has been done yet, make space
    if (!x->convolution)
        x->convolution = calloc(x->data_size, sizeof(float));

    // Inverse DFT
    // TODO: Potential optimization of saving this plan in the converb~ object
    fftwf_plan p2 = fftwf_plan_dft_c2r_1d(x->data_size, x->input_dft, x->convolution, FFTW_ESTIMATE);
    fftwf_execute(p2);
    fftwf_destroy_plan(p2);

    normalize(x->convolution, x->data_size);
}

void normalize(float* x, int size) {
    float max = 0;

    for (int i = 0; i < size; i++) {
        if (fabs(x[i]) > max)
            max = fabs(x[i]);
    }

    if (max != 0.0) {
        for (int i = 0; i < size; i++)
            x[i] /= max;
    }
}

void complex_multiply(float c1[2], float c2[2]) {
    c1[0] = (c1[0] * c2[0]) - (c1[1] * c2[1]);
    c1[1] = (c1[0] * c2[1]) + (c1[1] * c2[0]);
}