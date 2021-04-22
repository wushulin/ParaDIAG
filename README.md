# ParaDiag
ParaDiag is a class of Parallel-in-Time (PinT) algorithms based on the diagonalization technique. 
ParaDiag can handle both the dissipative problems (e.g., the advection-diffusion equations),
the purely wave equations (e.g., the acoustic equation), and optimal control of wave equations.
If you used the same technique in other applications, we would be happy to add your codes to this list.

# List of codes

- ParaDIAG_V1_for_ADE.m : the advection-diffusion equation by ParaDiag-I 
- ParaDIAG_V2_WR_for_ADE.m : the advection-diffusion equation by ParaDiag-II-Waveform Relaxation (WR) Variant
- ParaDIAG_V2_Parareal_for_ADE.m: the advection-diffusion equation by ParaDiag-II-Parareal Variant
- ParaDIAG_V2_GMRES_LinearWave_2D.m : the 2D wave equation by ParaDiag-II preconditioned GMRES
- ParaDIAG_V2_GMRES_LinearWaveOPT_2D.m: the optimal control of the wave equation by ParaDiag-II preconditioned GMRES
- ParaDIAG_V1_Hybrid_for_ADE.m: this the MATLAB code used for testing the new direct PinT algorithm (a newly updated material in Section 2.1.2)

- Parallel Codes.zip: parallel Fortran codes for the ParaDiag-II-Waveform Relaxation (WR) Variant
- ParaDIAG_V1_Hybrid_for_Wave.c: parallel codes based on C language for testing the new direct PinT algorithm (a newly updated material in Section 2.1.2)

# Citing reference (Please cite our work if appropriate)
- Martin J. Gander,  Jun Liu , Shu-Lin Wu, Xiaoqiang Yue, and Tao Zhou. 
"[ParaDiag: parallel-in-time algorithms based on the diagonalization technique.](https://arxiv.org/abs/2005.09158)" arXiv preprint arXiv:2005.09158 (2021).

# Current developers
- [Dr. Jun Liu](https://junliu2050.github.io/)
- [Dr. Shu-Lin Wu](https://www.researchgate.net/profile/Shu_Lin_Wu)
- [Dr. Xiaoqiang Yue](https://scholar.google.com/citations?user=oMMBhwgAAAAJ&hl=en)

# Other contributors
- [Dr. Martin J. Gander](https://www.unige.ch/~gander/)
- [Dr. Tao Zhou](http://lsec.cc.ac.cn/~tzhou/)

# Licence Information: 

This library is free software. 
It can be redistributed and/or modified under the terms of the [MIT License](https://opensource.org/licenses/MIT).

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Copyright (c) 2020 by Jun Liu, Shu-Lin Wu, and Xiaoqiang Yue.
