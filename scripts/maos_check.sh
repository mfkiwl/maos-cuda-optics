#!/bin/bash

#Sanity check maos
#NFIRAOS default
if [ -n "$1" ];then
    D=$1
    shift
    args="$@"
else
    D=30
    args=
fi
args+=" aper.d=$D"
case $D in
    2)
	#4/4/2017
	#REF=(90.1682 91.9371 166.001 133.074 133.682 140.137 130.599 138.509 131.078 589.324 225.003 588.795 220.674 117.774 120.16 588.795 211.419 150.461 135.644)
	#10/29/2018
	REF=(88.44 88.96 143.03 140.22 140.82 155.58 137.95 155.13 137.85 589.32 226.63 588.79 218.41 117.42 119.60 588.79 209.86 149.96 138.61)
	;;
    5)
	#6/7/2017
	REF=(107.10 108.90 137.96 132.96 132.58 132.10 132.73 132.85 121.69 510.76 259.97 104.77 107.74 125.48 130.55 121.69 122.34 253.28 136.87)
	;;
    10)
	#12/13/2018
	#REF=(108.80 109.09 127.83 129.54 130.46 129.81 error 129.03 115.42 353.55 287.17 104.30 104.76 139.21 143.33 195.46 165.91 217.42 128.84)
	#2/11/2019
	REF=(109.16 109.45 132.43 129.49 129.76 129.59 129.73 129.28 114.41 365.81 277.64 106.79 107.76 143.38 145.59 193.78 166.22 215.75 126.69)
	;;
    30)
	#9/14/2018
	#REF=(112.55 113.02 133.92 139.07 132.86 error 116.21 365.25 364.84 109.83 109.82 145.30 146.51 361.19 361.10 154.76 125.21)
	#12/13/2018: with sim.mffocus=1 as default
	#REF=(112.68 112.93 134.55 141.49 134.51 error 118.06 375.25 370.70 109.01 109.02 146.39 147.32 343.20 368.34 155.83 143.43)
	#2/11/2019
	REF=(112.68 112.93 134.58 141.23 133.94 146.73 117.34 375.55 370.68 108.96 109.09 146.35 147.25 343.93 367.85 156.07 138.63)

esac
fnlog=maos_check_${D}.log
fnres=maos_check_${D}.res
echo > $fnlog
echo "D is ${D}m. DM order is $((D*2))." |tee $fnres
ii=0

function run_maos(){
    kind="$1"
    shift
    echo "$kind" >>$fnlog
    stdbuf -o0 printf "%-20s" "$kind" | tee $fnres
    stdbuf -o0 ./maos sim.end=1000 $args "$*" >> $fnlog
    if [ $? == 0 ];then
	RMS[ii]=$(tail -n5 $fnlog |grep 'Mean:' |cut -d ':' -f 2)
	a=${RMS[$ii]%.*}
	ans=0
    else
	RMS[ii]='error'
	a=0
	ans=1
    fi

    b=${REF[$ii]%.*}
    stdbuf -o0 printf "%7s nm (%7s nm)" "${RMS[$ii]}" "${REF[$ii]}" |tee $fnres
    if [ $a -ne 0 -a "$b" != "error" ];then
	printf "%4s%%\n" $(((a-b)*100/b)) |tee $fnres
    else
	echo |tee $fnres
    fi
    ii=$((ii+1)) 
}



run_maos "Ideal fit (cpu): " sim.idealfit=1 

run_maos "Ideal tomo (cpu):" sim.idealtomo=1  

run_maos "LGS MCAO (inte): " recon.split=0 tomo.precond=0

run_maos "LGS MCAO (CG):   " tomo.precond=0

run_maos "LGS MCAO (FDPCG):" tomo.precond=1

run_maos "LGS MCAO (CBS):  " tomo.alg=0 fit.alg=0

if [ $D -le 10 ];then
run_maos "LGS MCAO (SVD):  " tomo.alg=2 fit.alg=2

run_maos "LGS MCAO (MVM):  " atmr.os=[2] tomo.precond=1 tomo.maxit=100 fit.alg=0 recon.mvm=1
fi

run_maos "LGS MOAO:        " evl.moao=0 moao.dx=[1/2]

run_maos "LGS GLAO (inte): " dm_single.conf  recon.glao=1 recon.split=0 wfs_lgs_ttf.conf

run_maos "LGS GLAO (split):" dm_single.conf  recon.glao=1 recon.split=1 wfs_lgs_ttf.conf

run_maos "NGS SCAO (inte): " -cscao_ngs.conf recon.split=0

run_maos "NGS SCAO (split):" -cscao_ngs.conf recon.split=1

run_maos "NGS MCAO (inte): " -cmcao_ngs.conf recon.split=0

run_maos "NGS MCAO (split):" -cmcao_ngs.conf recon.split=1

run_maos "SCAO LGS (inte): " -cscao_lgs.conf recon.split=0 powfs.dsa=[0.5 -1] 

run_maos "SCAO LGS (split):" -cscao_lgs.conf recon.split=1

run_maos "LGS LTAO (inte): " dm_single.conf fov_oa.conf recon.split=0 

run_maos "LGS LTAO (split):" dm_single.conf fov_oa.conf recon.split=1

echo ${RMS[*]} |tee $fnres

exit $ans
