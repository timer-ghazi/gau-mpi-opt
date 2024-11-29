
# Embarrassingly Parallel Gaussian Geometry Optimization with an External Code: High-Level Energies + Gaussian Optimizer
###  Q.T., Marquette U., 2018--2024

 * Currently supported electronic structure theory (EST) codes:
    - MRCC
    - Molpro

 * How to use:
   0. Make sure to copy this script into your work directory:
      cp /mmfs1/home/timerghazinq/local/gau-opt/G16-MPI-Grad-Opt-prep.sh .
      The script may require minimal customization such as molecular charge and
      multiplicity on top of the file.
   1. Run this script with a XYZ file name: ./G16-MPI-Grad-Opt-prep.sh mymolecule.xyz
      This will rotate the coordinates into standard orientation and also
      generate mymolecule.symm file with proposed symmetry constraints
   2. Correct any issues with generated constraints in mymolecule.symm
      E.g. any coordinate equal to 0.0 is assumed to be in a symmetry plane
      with zero force ('X00 = ZERO'). The algorithm is primitive--human editing is required!
      Make sure to remove any constraints for C₁ molecules!
   3. Run this script again, now with two files:
      ./G16-MPI-Grad-Opt-prep.sh mymolecule.xyz mymolecule.symm optional_job_name
      This will generate necessary input files and scripts setup for your system.
      If specified, optional_job_name will be used as the main file name.
       - To daisy-chain several optimizations (e.g. with increasing level of theory) use:
         ./G16-MPI-Grad-Opt-prep.sh mymolecule.xyz mymolecule.symm job1 job2 job3 ...
         This will create job1.com, job2.com, etc. files sequentially run from the
         mymolecule.sh script, each job starting with the previous job's checkpoint file.
   4. Run preliminary DFT frequency job dft-mymolecule-Freq.com as a regular Gaussian job
      This will provide a good starting Hessian for the geometry optimization.
   5. Edit mymolecule.com to setup the MOLPRO calculation and Gaussian optimization
      The MOLPRO input header is embedded in the Gaussian .com file as a comment, so
      each line starts with '!' that will be automatically removed later.
      Also, the MOLPRO input header is defined on top of this script (no '!' there).
   6. Adjust, if needed, and submit the Slurm job script mymolecule.sh
      - if at all possible, run only one MRCC job per node ('#SBATCH --ntasks-per-node=1')
      - monitor the messages in a Slurm output file mymolecule-*.SLURM for job progress/errors

  Note:
  * Only simple symmetry constraints implemented, but they can significantly cut the number of
    individual MOLPRO energy evaluations; automatically generated constraints will have mistakes,
    so editing the .symm file is always necessary.
    For C₁ molecules, make sure that there are no constraints (i.e. 'X01 = X01', etc. for all)

  * It is recommended to start with a lower coupled-cluster level optimization, and use the optimized
    geometry (and force constants) as a starting point for optimization with the next, higher level,
    which in turn will be the starting point for the higher yet theory-level optimization.
    E.g., to perform CCSDT(Q) optimization, a CCSD-->CCSD(T)-->CCSDT(Q) sequence can be used.
    Two options to do this:
     - After an optimization is complete, use copy2next-Gau-MRCC.sh script to produce the
       input & submittion script for the next step (.chk file from the previous job will be used)
     - Generate and run the entire sequence at once (see above):
       ./G16-MPI-Grad-Opt-prep.sh mymolecule.xyz mymolecule.symm mymolecule-CCSD mymolecule-CCSDt mymolecule-CCSDTq

  Requires:
   * Rest of this package--Python scripts in $PYSCRDIR
   * MRCC (in $MRCCDIR) and Gaussian (as a module)
   * Anaconda in $CONDADIR with the following packages installed:
      - Python 2.7 (installed as 'py2' environment)
      - numpy
      - mpi4py (which also installs mpich MPI library)
