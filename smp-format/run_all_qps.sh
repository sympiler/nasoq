#!/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --export=ALL
#SBATCH --job-name="spa_maros"
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes
#SBATCH --mail-user=kazem.cheshmi@gmail.com
#SBATCH --nodes=1
#SBATCH --output="all_qps.%j.%N.out"
#SBATCH --time 23:55:00

BINLIB=$1
PATHMAIN=$2
mode=$3
sol_mod=2
if [ "$#" -ge 4 ]; then
sol_mod=$4
fi

out_iter=3
in_iter=3
reg_diag=9
if [ "$#" -ge 5 ]; then
in_iter=$5
reg_diag=$6
fi

tol=15
eps=6

#1:only bbh;
#2: only maros;
#3: only contact;
#4: only cloth;
#10: all

#source /home/kazem/intel/mkl/bin/mklvars.sh intel64
#module load intel/2018.3
#export OMP_NUM_THREADS=20

echo "Tool Name,Problem Name,Hessian dim,Hessian NNZ,# of Eq Constraints,Eq Constraint NNZ,# of Ineq Const,Ineq Constraint NNZ,# of Threads,eps_abs,Outer GMRES Iter,Inner GMRES Iter,GMRES Tol,Diagonal Pert,Status,# of Iterations,Time (s),Active-set Size,Constraint Satisfaction Inf,Residual Lagrangian inf,Primal Obj,Dual Obj,Obj Value,Non-negativity Inf,Complementarity Inf,Problem Type,";



#for in_iter in {20,5,10}; do
    #for eps in {3,6,9}; do
#eps=3;
	    out_iter=$in_iter;
		#for reg_diag in {7..11}; do
		#for reg_diag in {9,10}; do
		#for reg_diag in {11,12}; do

########################## BBH
if [ $mode -eq 1 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/biharmonic_qps/
    ## bound harmonic
    #out_iter=2
    #in_iter=2
    #eps=6

    #for n in {hand,lion}; do
    n=hand
        for f in {1..3}; do
            $BINLIB 2 $PATHQP/${n}${f}_s_P.mtx $PATHQP/${n}${f}_s_q.mtx $PATHQP/${n}${f}_s_l.mtx $PATHQP/${n}${f}_s_A.mtx $PATHQP/${n}${f}_s_u.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod ${n} "libigl" "Image deformation"
            printf "Image deformation,"
            echo""
        done
    #done

    n=lion
        for f in {1..3}; do
            $BINLIB 2 $PATHQP/${n}${f}_s_P.mtx $PATHQP/${n}${f}_s_q.mtx $PATHQP/${n}${f}_s_l.mtx $PATHQP/${n}${f}_s_A.mtx $PATHQP/${n}${f}_s_u.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod ${n} "libigl" "Image deformation"
            printf "Image deformation,"
            echo""
        done

     for n in {alligator,woody}; do
        for f in {1..5}; do
            $BINLIB 2 $PATHQP/${n}${f}_P.mtx $PATHQP/${n}${f}_q.mtx $PATHQP/${n}${f}_l.mtx $PATHQP/${n}${f}_A.mtx $PATHQP/${n}${f}_u.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod ${n} "libigl" "Image deformation"
            printf "Image deformation,"
            echo""
        done
    done
	n=alligator
	for f in {6..9}; do
            $BINLIB 2 $PATHQP/${n}${f}_P.mtx $PATHQP/${n}${f}_q.mtx $PATHQP/${n}${f}_l.mtx $PATHQP/${n}${f}_A.mtx $PATHQP/${n}${f}_u.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod ${n} "libigl" "Image deformation"
            printf "Image deformation,"
            echo""
        done


    #for f in {10_15584,10_166,10_3900,10_643,15_15570,15_173,15_3917,15_629,20_15601,20_170,20_3921,20_640,2_15654,2_158,2_3871,2_616,3_15585,3_161,3_3912,3_618,4_155,4_15630,4_3859,4_622,5_15573,5_167,5_3938,5_618}; do
    #    echo $BINLIB 2 $PATHQP/bbh_square_matrix_${f}.mtx $PATHQP/bbh_square_bounds_${f}_q.txt $PATHQP/bbh_square_bounds_${f}_l.txt $PATHQP/bbh_square_bounds_${f}_A.mtx $PATHQP/bbh_square_bounds_${f}_u.txt $reg_diag $out_iter $in_iter $eps $tol
    #    printf "Image deformation,"
    #    echo""
    #done
    echo ""
fi



######################### Maros Mezaros
if [ $mode -eq 2 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/MarosMezaros/
 #{AUG2D ,AUG2DC ,AUG2DCQP ,AUG2DQP ,AUG3D ,AUG3DC ,AUG3DCQP ,AUG3DQP ,BOYD1 ,BOYD2 ,CONT-050 ,CONT-100 ,CONT-101 ,CONT-200 ,CONT-201 ,CONT-300 ,CVXQP1_L ,CVXQP1_M ,CVXQP1_S ,CVXQP2_L ,CVXQP2_M ,CVXQP2_S ,CVXQP3_L ,CVXQP3_M ,CVXQP3_S ,DPKLO1 ,DTOC3 ,DUAL1 ,DUAL2 ,DUAL3 ,DUAL4 ,DUALC1 ,DUALC2 ,DUALC5 ,DUALC8 ,EXDATA ,GENHS28 ,GOULDQP2 ,GOULDQP3 ,HS118 ,HS21 ,HS268 ,HS35 ,HS35MOD ,HS51 ,HS52 ,HS53 ,HS76 ,HUES-MOD ,HUESTIS ,KSIP ,LASER ,LISWET1 ,LISWET10 ,LISWET11 ,LISWET12 ,LISWET2 ,LISWET3 ,LISWET4 ,LISWET5 ,LISWET6 ,LISWET7 ,LISWET8 ,LISWET9 ,LOTSCHD ,MOSARQP1 ,MOSARQP2 ,POWELL20 ,PRIMAL1 ,PRIMAL2 ,PRIMAL3 ,PRIMAL4 ,PRIMALC1 ,PRIMALC2 ,PRIMALC5 ,PRIMALC8 ,Q25FV47 ,QADLITTL ,QAFIRO ,QBANDM ,QBEACONF ,QBORE3D ,QBRANDY ,QCAPRI ,QE226 ,QETAMACR ,QFFFFF80 ,QFORPLAN ,QGFRDXPN ,QGROW15 ,QGROW22 ,QGROW7 ,QISRAEL ,QPCBLEND ,QPCBOEI1 ,QPCBOEI2 ,QPCSTAIR ,QPILOTNO ,QPTEST ,QRECIPE ,QSC205 ,QSCAGR25 ,QSCAGR7 ,QSCFXM1 ,QSCFXM2 ,QSCFXM3 ,QSCORPIO ,QSCRS8 ,QSCSD1 ,QSCSD6 ,QSCSD8 ,QSCTAP1 ,QSCTAP2 ,QSCTAP3 ,QSEBA ,QSHARE1B ,QSHARE2B ,QSHELL ,QSHIP04L ,QSHIP04S ,QSHIP08L ,QSHIP08S ,QSHIP12L ,QSHIP12S ,QSIERRA ,QSTAIR ,QSTANDAT ,S268 ,STADAT1 ,STADAT2 ,STADAT3 ,STCQP1 ,STCQP2 ,TAME ,UBH1 ,VALUES ,YAO ,ZECEVIC2 }

 ## convex
### 3D problems
#for f in {AUG2D,AUG2DC,AUG2DCQP,AUG2DQP,AUG3D,AUG3DC,AUG3DCQP,AUG3DQP,BOYD1,BOYD2,CONT-050,CONT-100,CONT-101,CONT-200,CONT-201,CONT-300,CVXQP1_L,CVXQP1_M,CVXQP1_S,CVXQP2_L,CVXQP2_M,CVXQP2_S,CVXQP3_L,CVXQP3_M,CVXQP3_S,DPKLO1,DTOC3,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC2,DUALC5,DUALC8,EXDATA,GENHS28,GOULDQP2,GOULDQP3,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,LOTSCHD,MOSARQP1,MOSARQP2,POWELL20,PRIMAL1,PRIMAL2,PRIMAL3,PRIMAL4,PRIMALC1,PRIMALC2,PRIMALC5,PRIMALC8,Q25FV47,QADLITTL,QAFIRO,QBANDM,QBEACONF,QBORE3D,QBRANDY,QCAPRI,QE226,QETAMACR,QFFFFF80,QFORPLAN,QGFRDXPN,QGROW15,QGROW22,QGROW7,QISRAEL,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPILOTNO,QPTEST,QRECIPE,QSC205,QSCAGR25,QSCAGR7,QSCFXM1,QSCFXM2,QSCFXM3,QSCORPIO,QSCRS8,QSCSD1,QSCSD6,QSCSD8,QSCTAP1,QSCTAP2,QSCTAP3,QSEBA,QSHARE1B,QSHARE2B,QSHELL,QSHIP04L,QSHIP04S,QSHIP08L,QSHIP08S,QSHIP12L,QSHIP12S,QSIERRA,QSTAIR,QSTANDAT,README,S268,STADAT1,STADAT2,STADAT3,STCQP1,STCQP2,TAME,UBH1,VALUES,YAO,ZECEVIC2,cleanup_dir,extract_sif,sif2mat,AUG2D,AUG2DC,AUG2DCQP,AUG2DQP,AUG3D,AUG3DC,AUG3DCQP,AUG3DQP,BOYD1,BOYD2,CONT-050,CONT-100,CONT-101,CONT-200,CONT-201,CONT-300,CVXQP1_L,CVXQP1_M,CVXQP1_S,CVXQP2_L,CVXQP2_M,CVXQP2_S,CVXQP3_L,CVXQP3_M,CVXQP3_S,DPKLO1,DTOC3,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC2,DUALC5,DUALC8,EXDATA,GENHS28,GOULDQP2,GOULDQP3,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,LOTSCHD,MOSARQP1,MOSARQP2,POWELL20,PRIMAL1,PRIMAL2,PRIMAL3,PRIMAL4,PRIMALC1,PRIMALC2,PRIMALC5,PRIMALC8,Q25FV47,QADLITTL,QAFIRO,QBANDM,QBEACONF,QBORE3D,QBRANDY,QCAPRI,QE226,QETAMACR,QFFFFF80,QFORPLAN,QGFRDXPN,QGROW15,QGROW22,QGROW7,QISRAEL,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPILOTNO,QPTEST,QRECIPE,QSC205,QSCAGR25,QSCAGR7,QSCFXM1,QSCFXM2,QSCFXM3,QSCORPIO,QSCRS8,QSCSD1,QSCSD6,QSCSD8,QSCTAP1,QSCTAP2,QSCTAP3,QSEBA,QSHARE1B,QSHARE2B,QSHELL,QSHIP04L,QSHIP04S,QSHIP08L,QSHIP08S,QSHIP12L,QSHIP12S,QSIERRA,QSTAIR,QSTANDAT,S268,STADAT1,STADAT2,STADAT3,STCQP1,STCQP2,TAME,UBH1,VALUES,YAO,ZECEVIC2};
#for f in {AUG2D,AUG2DC,AUG2DCQP,AUG2DQP,AUG3D,AUG3DC,AUG3DCQP,AUG3DQP}; #requires 1 iter refinement except for last
#for f in {AUG2DC,AUG2DCQP,AUG3DC,AUG3DCQP,CONT-050,CONT-100,CONT-200,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC5,HS118,HS21,HS268,HS35,HS35MOD,HS52,HS76,HUES-MOD,HUESTIS,KSIP,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,MOSARQP1,MOSARQP2,POWELL20,PRIMAL3,PRIMAL4,PRIMALC1,PRIMALC8,QGFRDXPN,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPTEST,QSCAGR25,QSCFXM3,QSEBA,S268,STADAT1,STCQP1,STCQP2,YAO};
#for f in {BOYD1,BOYD2};
#for f in {CVXQP1_L,CVXQP1_M,CVXQP1_S,CVXQP2_L,CVXQP2_M,CVXQP2_S,CVXQP3_L,CVXQP3_M,CVXQP3_S};
#for f in {DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC2,DUALC5,DUALC8};
#for f in {GENHS28,GOULDQP2,GOULDQP3};
#for f in {HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS};
#for f in {LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,LOTSCHD};

# all convex problems
#for out_iter in {1..15}; do
#	for in_iter in {2,4,6,8,10,15,20}; do
#		for reg_diag in {6,8,10,12,14}; do
#			for tol in {15,12,8,10,18}; do
#				for eps in {5,6,10}; do
#					for f in {AUG2DC,AUG2DCQP,AUG3DC,AUG3DCQP,BOYD1,CONT-050,CONT-100,CONT-200,CONT-201,CONT-300,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC5,GENHS28,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,MOSARQP1,MOSARQP2,POWELL20,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPTEST,S268,STCQP1,STCQP2,TAME,YAO};
##for f in {MOSARQP1,POWELL20,YAO,QPCSTAIR,MOSARQP2};
#					do
#					 $BINLIB $PATHQP/${f}_P $PATHQP/${f}_q $PATHQP/${f}_l $PATHQP/${f}_A $PATHQP/${f}_u $reg_diag $out_iter $in_iter $eps $tol;
#					 echo""
#					done
#				done
#			done
#		done
#	done
#done

#out_iter=3
#in_iter=3
#eps=6
#for f in {AUG2DC,AUG2DCQP,AUG3DC,AUG3DCQP,BOYD1,CONT-050,CONT-100,CONT-200,CONT-201,CONT-300,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC5,GENHS28,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,MOSARQP1,MOSARQP2,POWELL20,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPTEST,S268,STCQP1,STCQP2,TAME,YAO};
#for f in {AUG2D,AUG2DC,AUG2DCQP,AUG2DQP,AUG3D,AUG3DC,AUG3DCQP,AUG3DQP,BOYD1,BOYD2,CONT-050,CONT-100,CONT-101,CONT-200,CONT-201,CONT-300,CVXQP1_L,CVXQP1_M,CVXQP1_S,CVXQP2_L,CVXQP2_M,CVXQP2_S,CVXQP3_L,CVXQP3_M,CVXQP3_S,DPKLO1,DTOC3,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC2,DUALC5,DUALC8,EXDATA,GENHS28,GOULDQP2,GOULDQP3,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,LOTSCHD,MOSARQP1,MOSARQP2,POWELL20,PRIMAL1,PRIMAL2,PRIMAL3,PRIMAL4,PRIMALC1,PRIMALC2,PRIMALC5,PRIMALC8,Q25FV47,QADLITTL,QAFIRO,QBANDM,QBEACONF,QBORE3D,QBRANDY,QCAPRI,QE226,QETAMACR,QFFFFF80,QFORPLAN,QGFRDXPN,QGROW15,QGROW22,QGROW7,QISRAEL,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPILOTNO,QPTEST,QRECIPE,QSC205,QSCAGR25,QSCAGR7,QSCFXM1,QSCFXM2,QSCFXM3,QSCORPIO,QSCRS8,QSCSD1,QSCSD6,QSCSD8,QSCTAP1,QSCTAP2,QSCTAP3,QSEBA,QSHARE1B,QSHARE2B,QSHELL,QSHIP04L,QSHIP04S,QSHIP08L,QSHIP08S,QSHIP12L,QSHIP12S,QSIERRA,QSTAIR,QSTANDAT,S268,STADAT1,STADAT2,STADAT3,STCQP1,STCQP2,TAME,UBH1,VALUES,YAO,ZECEVIC2};
for f in {AUG2DC,AUG2DCQP,AUG3DC,AUG3DCQP,CONT-050,CONT-100,CONT-200,CONT-201,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC5,GENHS28,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,MOSARQP1,MOSARQP2,POWELL20,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPTEST,S268,STCQP1,STCQP2,TAME,YAO};
#for f in {LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,YAO};
do
  $BINLIB 2 $PATHQP/${f}_P $PATHQP/${f}_q $PATHQP/${f}_l $PATHQP/${f}_A $PATHQP/${f}_u $reg_diag $out_iter $in_iter $eps $tol $sol_mod "MarosMeszaros" "Maros-Meszaros" "Maros-Meszaros";
 printf "Maros Mezaros,"
 echo""
done

for f in {LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9};
do
  $BINLIB 2 $PATHQP/${f}_P $PATHQP/${f}_q $PATHQP/${f}_l $PATHQP/${f}_A $PATHQP/${f}_u $reg_diag $out_iter $in_iter $eps $tol $sol_mod "MarosMeszaros" "Maros-Meszaros" "Maros-Meszaros";
 printf "Maros Mezaros,"
 echo""
done
#echo""
#out_iter=5
#in_iter=5
#eps=6
##for f in {AUG2DC,AUG2DCQP,AUG3DC,AUG3DCQP,BOYD1,CONT-050,CONT-100,CONT-200,CONT-201,CONT-300,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC5,GENHS28,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,MOSARQP1,MOSARQP2,POWELL20,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPTEST,S268,STCQP1,STCQP2,TAME,YAO};
#for f in {AUG2DC,AUG2DCQP,AUG3DC,AUG3DCQP,CONT-050,CONT-100,CONT-200,CONT-201,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC5,GENHS28,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,MOSARQP1,MOSARQP2,POWELL20,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPTEST,S268,STCQP1,STCQP2,TAME,YAO};
##for f in {LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,YAO};
#do
# $BINLIB 2 $PATHQP/${f}_P $PATHQP/${f}_q $PATHQP/${f}_l $PATHQP/${f}_A $PATHQP/${f}_u $reg_diag $out_iter $in_iter $eps $tol;
# echo""
#done

#for f in {AUG2DC,AUG2DCQP,AUG3DC,AUG3DCQP,BOYD1,CONT-050,CONT-100,CONT-200,CONT-201,CONT-300,DUAL1,DUAL2,DUAL3,DUAL4,DUALC1,DUALC5,GENHS28,HS118,HS21,HS268,HS35,HS35MOD,HS51,HS52,HS53,HS76,HUES-MOD,HUESTIS,KSIP,LASER,LISWET1,LISWET10,LISWET11,LISWET12,LISWET2,LISWET3,LISWET4,LISWET5,LISWET6,LISWET7,LISWET8,LISWET9,MOSARQP1,MOSARQP2,POWELL20,QPCBLEND,QPCBOEI1,QPCBOEI2,QPCSTAIR,QPTEST,S268,STCQP1,STCQP2,TAME,YAO}
echo""
fi


######################### Contact problems
if [ $mode -eq 3 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/MarosMezaros/
# "2D Problems"
$BINLIB 1 $PATHQP/multiple_contacts_hP.mtx $PATHQP/multiple_contacts_q none none $PATHQP/multiple_contacts_A.mtx $PATHQP/multiple_contacts_l $reg_diag $out_iter $in_iter $eps $tol $sol_mod  "Simulation2D" "NASOQ" "Contact Simulation" ;
printf "Contact simulation,"
echo ""


$BINLIB 1 $PATHQP/contact_starting_hP.mtx $PATHQP/contact_starting_q none none $PATHQP/contact_starting_A.mtx $PATHQP/contact_starting_l $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Simulation2D" "NASOQ" "Contact Simulation" ;
printf "Contact simulation,"
echo ""

$BINLIB 1 $PATHQP/no_contact_hP.mtx $PATHQP/no_contact_q none none $PATHQP/no_contact_A.mtx $PATHQP/no_contact_l $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Simulation2D" "NASOQ" "Contact Simulation" ;
printf "Contact simulation,"
echo ""

$BINLIB 1 $PATHQP/recontactqps_hP.mtx $PATHQP/recontactqps_q none none $PATHQP/recontactqps_A.mtx $PATHQP/recontactqps_l $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Simulation2D" "NASOQ" "Contact Simulation" ;
printf "Contact simulation,"
echo""


### 3D problems
for f in {8..10};
do
  ${BINLIB} 1 $PATHQP/Ph${f}.mtx $PATHQP/q${f} none none $PATHQP/A${f} $PATHQP/l${f} $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Simulation3D" "NASOQ" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

#for f in {176..201};
for f in {176..178};
do
 ${BINLIB} 1 $PATHQP/Ph${f}.mtx $PATHQP/q${f} none none $PATHQP/A${f} $PATHQP/l${f} $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Simulation3D" "NASOQ" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

echo ""
fi


######################### cloth problems
if [ $mode -eq 4 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/cloth_sim/
#out_iter=1
#in_iter=1
#box
#for f in {1..99};
#do
# ${BINLIB} 1 $PATHQP/box_${f}_P.mtx $PATHQP/box_${f}_q.mtx $PATHQP/box_${f}_A.mtx $PATHQP/box_${f}_b.mtx $PATHQP/box_${f}_C.mtx $PATHQP/box_${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# echo ""
#done

#boxmove
#for f in {1..16};
#do
# ${BINLIB} 1 $PATHQP/box_move_${f}_P.mtx $PATHQP/box_move_${f}_q.mtx $PATHQP/box_move_${f}_A.mtx $PATHQP/box_move_${f}_b.mtx $PATHQP/box_move_${f}_C.mtx $PATHQP/box_move_${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# echo ""
#done

#boxstill
#for f in {1..17};
#do
# ${BINLIB} 1 $PATHQP/box_still_${f}_P.mtx $PATHQP/box_still_${f}_q.mtx $PATHQP/box_still_${f}_A.mtx $PATHQP/box_still_${f}_b.mtx $PATHQP/box_still_${f}_C.mtx $PATHQP/box_still_${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# echo ""
#done

# nail
for f in {7..9};
do
 ${BINLIB} 1 $PATHQP/nail0${f}t_hP.mtx $PATHQP/nail0${f}t_q.mtx $PATHQP/nail0${f}t_A.mtx $PATHQP/nail0${f}t_b.mtx $PATHQP/nail0${f}t_C.mtx $PATHQP/nail0${f}t_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Nail" "eol-cloth repository" "Cloth Simulation";
  printf "Cloth simulation,"
  echo ""
done

#for f in {10..18};
#do
# ${BINLIB} 1 $PATHQP/nail${f}t_hP.mtx $PATHQP/nail${f}t_q.mtx $PATHQP/nail${f}t_A.mtx $PATHQP/nail${f}t_b.mtx $PATHQP/nail${f}t_C.mtx $PATHQP/nail${f}t_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Nail" "eol-cloth repository" "Cloth Simulation";;
#  printf "Cloth simulation,"
#  echo ""
#done

#wind
#for f in {1..74};
#do
# ${BINLIB} 1 $PATHQP/wind_${f}_P.mtx $PATHQP/wind_${f}_q.mtx $PATHQP/wind_${f}_A.mtx $PATHQP/wind_${f}_b.mtx $PATHQP/wind_${f}_C.mtx $PATHQP/wind_${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# echo ""
#done

echo ""
fi

######################### Contact problems
if [ $mode -eq 5 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/gauss_qps/
for f in {0..1};
do
 ${BINLIB} 1 $PATHQP/bunny2_bunny${f}_hP.mtx $PATHQP/bunny2_bunny${f}_q.mtx none none $PATHQP/bunny2_bunny${f}_C.mtx $PATHQP/bunny2_bunny${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Bunny" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done
#for f in {0..1};
#do
f=0
 ${BINLIB} 1 $PATHQP/bunny2_falling${f}_hP.mtx $PATHQP/bunny2_falling${f}_q.mtx none none $PATHQP/bunny2_falling${f}_C.mtx $PATHQP/bunny2_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Bunny" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
#done
for f in {0..1};
do
 ${BINLIB} 1 $PATHQP/bunny_bunny${f}_hP.mtx $PATHQP/bunny_bunny${f}_q.mtx none none $PATHQP/bunny_bunny${f}_C.mtx $PATHQP/bunny_bunny${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Bunny" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done
#for f in {0..1};
#do
f=0
 ${BINLIB} 1 $PATHQP/bunny_falling${f}_hP.mtx $PATHQP/bunny_falling${f}_q.mtx none none $PATHQP/bunny_falling${f}_C.mtx $PATHQP/bunny_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Bunny" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
#done

for f in {0..4};
do
 ${BINLIB} 1 $PATHQP/horse_horse${f}_hP.mtx $PATHQP/horse_horse${f}_q.mtx none none $PATHQP/horse_horse${f}_C.mtx $PATHQP/horse_horse${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Horse" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

#for f in {0..1};
#do
f=0
 ${BINLIB} 1 $PATHQP/horse_falling${f}_hP.mtx $PATHQP/horse_falling${f}_q.mtx none none $PATHQP/horse_falling${f}_C.mtx $PATHQP/horse_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Horse" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
#done

for f in {0..19};
do
 ${BINLIB} 1 $PATHQP/knight_knight${f}_hP.mtx $PATHQP/knight_knight${f}_q.mtx none none $PATHQP/knight_knight${f}_C.mtx $PATHQP/knight_knight${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Knight" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done
for f in {0..1};
do
 ${BINLIB} 1 $PATHQP/knight_falling${f}_hP.mtx $PATHQP/knight_falling${f}_q.mtx none none $PATHQP/knight_falling${f}_C.mtx $PATHQP/knight_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Knight" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

for f in {0..36};
do
 ${BINLIB} 1 $PATHQP/lamb_lamb${f}_hP.mtx $PATHQP/lamb_lamb${f}_q.mtx none none $PATHQP/lamb_lamb${f}_C.mtx $PATHQP/lamb_lamb${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Lamb" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done
#for f in {0..1};
#do
f=0
 ${BINLIB} 1 $PATHQP/lamb_falling${f}_hP.mtx $PATHQP/lamb_falling${f}_q.mtx none none $PATHQP/lamb_falling${f}_C.mtx $PATHQP/lamb_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Lamb" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
#done

#${BINLIB} 1 $PATHQP/last_qp_hP.mtx $PATHQP/last_qp_q.mtx none none $PATHQP/last_qp_C.mtx $PATHQP/last_qp_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
#printf "Contact simulation,"
#echo ""

for f in {1..3};
do
 ${BINLIB} 1 $PATHQP/matt_qp${f}_hP.mtx $PATHQP/matt_qp${f}_q.mtx none none $PATHQP/matt_qp${f}_C.mtx $PATHQP/matt_qp${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Matt" "Contact simulator" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

#for f in {0..40};
#do
# ${BINLIB} 1 $PATHQP/screwdriver_falling${f}_hP.mtx $PATHQP/screwdriver_falling${f}_q.mtx none none $PATHQP/screwdriver_falling${f}_C.mtx $PATHQP/screwdriver_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# printf "Contact simulation,"
# echo ""
#done



for f in {0..5};
do
 ${BINLIB} 1 $PATHQP/wolf_falling${f}_hP.mtx $PATHQP/wolf_falling${f}_q.mtx none none $PATHQP/wolf_falling${f}_C.mtx $PATHQP/wolf_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Wolf" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done


${BINLIB} 1 $PATHQP/garg3_garg30_hP.mtx $PATHQP/garg3_garg30_q.mtx none none $PATHQP/garg3_garg30_C.mtx $PATHQP/garg3_garg30_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Gargoyle" "GAUSS Repository" "Contact Simulation";
printf "Contact simulation,"
echo ""

f=0
#for f in {0..1};
#do
 ${BINLIB} 1 $PATHQP/garg3_floor${f}_hP.mtx $PATHQP/garg3_floor${f}_q.mtx none none $PATHQP/garg3_floor${f}_C.mtx $PATHQP/garg3_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Gargoyle" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
#done

for f in {0..14};
do
 ${BINLIB} 1 $PATHQP/garg2_garg2${f}_hP.mtx $PATHQP/garg2_garg2${f}_q.mtx none none $PATHQP/garg2_garg2${f}_C.mtx $PATHQP/garg2_garg2${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Gargoyle" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

for f in {0..87};
do
 ${BINLIB} 1 $PATHQP/garg_floor${f}_hP.mtx $PATHQP/garg_floor${f}_q.mtx none none $PATHQP/garg_floor${f}_C.mtx $PATHQP/garg_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Gargoyle" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done


for f in {0..15};
do
  ${BINLIB} 1 $PATHQP/arma1_arma${f}_hP.mtx $PATHQP/arma1_arma${f}_q.mtx none none $PATHQP/arma1_arma${f}_C.mtx $PATHQP/arma1_arma${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Armadillo" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

#for f in {0..2};
#do
#  ${BINLIB} 1 $PATHQP/arma1_floor${f}_hP.mtx $PATHQP/arma1_floor${f}_q.mtx none none $PATHQP/arma1_floor${f}_C.mtx $PATHQP/arma1_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# echo ""
#done

#for f in {0..4};
#do
#  ${BINLIB} 1 $PATHQP/arma2_floor${f}_hP.mtx $PATHQP/arma2_floor${f}_q.mtx none none $PATHQP/arma2_floor${f}_C.mtx $PATHQP/arma2_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# echo ""
#done

for f in {2..4};
do
  ${BINLIB} 1 $PATHQP/Bar2k_floor${f}_hP.mtx $PATHQP/Bar2k_floor${f}_q.mtx none none $PATHQP/Bar2k_floor${f}_C.mtx $PATHQP/Bar2k_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Bar2k" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

f=0
${BINLIB} 1 $PATHQP/Beam_floor${f}_hP.mtx $PATHQP/Beam_floor${f}_q.mtx none none $PATHQP/Beam_floor${f}_C.mtx $PATHQP/Beam_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Beam" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""

for f in {130..289};
do
  ${BINLIB} 1 $PATHQP/Beam_floor${f}_hP.mtx $PATHQP/Beam_floor${f}_q.mtx none none $PATHQP/Beam_floor${f}_C.mtx $PATHQP/Beam_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Beam" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

f=338
  ${BINLIB} 1 $PATHQP/Beam_floor${f}_hP.mtx $PATHQP/Beam_floor${f}_q.mtx none none $PATHQP/Beam_floor${f}_C.mtx $PATHQP/Beam_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Beam" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""

#for f in {0..47};
#do
#  ${BINLIB} 1 $PATHQP/brick0_brick${f}_hP.mtx $PATHQP/brick0_brick${f}_q.mtx none none $PATHQP/brick0_brick${f}_C.mtx $PATHQP/brick0_brick${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# printf "Contact simulation,"
# echo ""
#done

#for f in {0..3};
#do
#  ${BINLIB} 1 $PATHQP/cactus_floor${f}_hP.mtx $PATHQP/cactus_floor${f}_q.mtx none none $PATHQP/cactus_floor${f}_C.mtx $PATHQP/cactus_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod;
# echo ""
#done


for f in {2..524};
do
  ${BINLIB} 1 $PATHQP/Cube_floor${f}_hP.mtx $PATHQP/Cube_floor${f}_q.mtx none none $PATHQP/Cube_floor${f}_C.mtx $PATHQP/Cube_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Cube" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done

for f in {0..71};
do
  ${BINLIB} 1 $PATHQP/gargNested2_falling${f}_hP.mtx $PATHQP/gargNested2_falling${f}_q.mtx none none $PATHQP/gargNested2_falling${f}_C.mtx $PATHQP/gargNested2_falling${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Gargoyle2" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done


for f in {0..279};
do
  ${BINLIB} 1 $PATHQP/gargNested2_floor${f}_hP.mtx $PATHQP/gargNested2_floor${f}_q.mtx none none $PATHQP/gargNested2_floor${f}_C.mtx $PATHQP/gargNested2_floor${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Gargoyle2" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done


for f in {0..56};
do
  ${BINLIB} 1 $PATHQP/gargNested2_gargNested2${f}_hP.mtx $PATHQP/gargNested2_gargNested2${f}_q.mtx none none $PATHQP/gargNested2_gargNested2${f}_C.mtx $PATHQP/gargNested2_gargNested2${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Gargoyle2" "GAUSS Repository" "Contact Simulation";
 printf "Contact simulation,"
 echo ""
done


echo ""
fi

######################### Distortion mapping
if [ $mode -eq 6 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/gauss_qps/
#for f in {0..3};
#do
# ${BINLIB} 1 $PATHQP/InjBnd_Mapping${f}_hP.mtx $PATHQP/InjBnd_Mapping${f}_q.mtx $PATHQP/InjBnd_Mapping${f}_A.mtx $PATHQP/InjBnd_Mapping${f}_b.mtx $PATHQP/InjBnd_Mapping${f}_C.mtx $PATHQP/InjBnd_Mapping${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Elephant" "Distortion paper" "Distortion Mapping";
#  printf "Mesh processing,"
#  echo ""
#done


for f in {0..2};
do
 ${BINLIB} 1 $PATHQP/InjBnd_MappingImp${f}_hP.mtx $PATHQP/InjBnd_MappingImp${f}_q.mtx $PATHQP/InjBnd_MappingImp${f}_A.mtx $PATHQP/InjBnd_MappingImp${f}_b.mtx $PATHQP/InjBnd_MappingImp${f}_C.mtx $PATHQP/InjBnd_MappingImp${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Elephant" "Distortion paper" "Distortion Mapping";
  printf "Mesh processing,"
  echo ""
done

echo ""
fi

######################### Model construction
if [ $mode -eq 7 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/model_reconstruction/
for f in "bunny" "elephant" "farmer" "lamb" "unicorn" "wolf";
do
  ${BINLIB} 1 $PATHQP/${f}_hP.mtx $PATHQP/${f}_q.mtx $PATHQP/${f}_A.mtx $PATHQP/${f}_b.mtx $PATHQP/${f}_C.mtx $PATHQP/${f}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "AddingDepth" "Daniel and Marek paper" "Model Construction";
  printf "Model construction,"
  echo ""
done

echo ""
fi


######################### MPC
if [ $mode -eq 8 ] || [ $mode -eq 120 ]; then
PATHQP=$PATHMAIN/mpc_qps/
for f in "test01" "test02" "test03" "test04" "test05" "test06";do
for j in {0..19}; do
  ${BINLIB} 1 $PATHQP/${f}_${j}_P.mtx $PATHQP/${f}_${j}_q.mtx $PATHQP/${f}_${j}_A.mtx $PATHQP/${f}_${j}_b.mtx $PATHQP/${f}_${j}_C.mtx $PATHQP/${f}_${j}_d.mtx $reg_diag $out_iter $in_iter $eps $tol $sol_mod "Mpclib" "MPC Library" "Model Predictive Control";
  printf "MPC,"
  echo ""
done
done


echo ""
fi


#done
#done
#done
