#! /bin/bash
source /home/chenjiyuan/conda.env

energy=400
particle="e-"

if [ $particle = "e-" ]
then
    job=$((energy / 10))
elif [ $particle = "pi-" ]
then
    job=$((energy / 10 + 200))
fi

filename="/lustre/collider/chenjiyuan/dss-pid/run/e-pi-/root/training/job${job}_${particle}_${energy}MeV/hit_sel_${particle}_${energy}MeV_ana.root"
#tree=dp
event_index=2024
save_dir="/lustre/collider/chenjiyuan/dss-pid/figs/"
output="EventDisplay_${particle}_${energy}MeV.pdf"
#show=0

title="${energy}-MeV "
if [ $particle = "e-" ]
then
    title+='$e^-$'
elif [ $particle = "mu-" ]
then
    title+='$\mu^-$'
elif [ $particle = "pi-" ]
then
    title+='$\pi^-$'
elif [ $particle = "e+" ]
then
    title+='$e^+$'
elif [ $particle = "mu+" ]
then
    title+='$\mu^+$'
elif [ $particle = "pi+" ]
then
    title+='$\pi^+$'
elif [ $particle = "gamma" ]
then
    title+='$\gamma$'
elif [ $particle = "neutron" ]
then
    title+='$n$'
elif [ $particle = "proton" ]
then
    title+='$p$'
fi

#python display.py -f=$filename -i="$title" -e=$event_index -d=$save_dir -o=$output -s=$show
python display.py -f=$filename -i="$title" -e=$event_index -d=$save_dir -o="$output"
