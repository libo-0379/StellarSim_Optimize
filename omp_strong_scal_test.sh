#!/bin/bash
rootPath=`pwd`
workZone=omp_strong_scal_test
exeBin=/home/user/exe
inputFile=/home/user/input_files

for nthreads in 1 2 4 6 12 24 32
do
    workDir=${workZone}/nthreads${nthreads}
    mkdir -p $workDir
    cp -rf $inputFile $workDir
    cd $workDir

cat > yhrun.sh << EOF
#!/bin/bash
export OMP_NUM_THREADS=$nthreads
yhrun -N 1 -n 1 -c $nthreads $exeBin
EOF

    yhbatch -p thcp1 -N 1 ./yhrun.sh

    cd $rootPath
done
