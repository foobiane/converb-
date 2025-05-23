## converb~: A basic convolution reverb external in Pd 

**This project is still a work-in-progress. Expect bugs, malfunctions, and generally bad code.**

*Created by Ian Doherty, January 2025 - February 2025*

This Pure Data external convolves an input source signal with a preloaded impulse response file in real-time. The convolution is *not partitioned*: the delay in the output signal is thus a direct result of the length of the impulse response. Fast convolution (point-wise multiplication of FFTs) is used to minimize computational cost.

This external was created as part of my senior capstone project at the University of Illinois.

Potential improvements for future iterations of this external:
* Use partitioning to reduce IR delay.
* Use a more sophisticated vector multiplication algorithm to speed up the convolution step. This could consist of SIMD and/or alignment tricks in C.
* Investigate if FFTW plans can be saved and reused in the object proper.

See below for installation instructions. For usage, refer to `converb~-help.pd`. 

## Installation
1. If you haven't already, install PureData [here](https://puredata.info/).
2. Clone this repo using `git clone` on a command line or download the source code manually.
3. Install FFTW3 [here](https://fftw.org/download.html) and modify the Makefile to suit your particular system.
   * Windows: Ensure that you have [mingw-64w](https://www.mingw-w64.org/) installed on your system. Download the 64-bit precompiled version of FFTW. Place `fftw3.h` in `msys64/usr/local/include` and `libfftw3f-3.dll` in `msys64/usr/local/lib`. The Makefile is already configured for Windows.
   * OSX/Linux: Download the FFTW tarball from their website. In the command line, navigate to the download folder and run `./configure CFLAGS="-fPIC" --enable-float` for OSX or `./configure CFLAGS="-arch i386 -arch x86_64" --enable-float` for Linux. Next, run `make` and `sudo make install`. Finally, comment out the Windows line in the Makefile and uncomment the MacOSX / Linux line.
5. In the installation folder of convolve~, run `make` on the command line (Windows users, use mingw-64w; see above). For a system-wide installation, run `make install`.
6. Open `converb~-help.pd` to learn the usage for converb~. If `make install` was used, you will be able to add converb~ objects into any PureData instance.
	* *To view this in PureData, go to Help > Browser. converb~ should be listed.*

## Credits & Libraries Used
* [PureData](https://puredata.info/): Everything!
* [FFTW3](https://github.com/FFTW/fftw3): DFT generation
* [convolve_tilde](https://github.com/wbrent/convolve_tilde): Makefile design, FFTW3 usage, installation instructions
* [pure-data/pd-lib-builder](https://github.com/pure-data/pd-lib-builder): Makefile generation
* [bed~ by Eric Lyon](https://www.amazon.com/Designing-Audio-Objects-Max-MSP/dp/B009LLXIVC): Handling garrays in the Pd external lib
* [pure-data/externals-howto](https://github.com/pure-data/externals-howto?tab=readme-ov-file#atom-string): Helpful guide used during the creation of this external
