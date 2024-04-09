#! /bin/bash
source /home/chenjiyuan/conda.env

#Available modes: "display" for event display, and "projection" for energy projection. They must correspond to the names of the existing python files!
mode="projection"
energy=2000
particle="e-"

if [ $particle = "e-" ]
then
    job=$((energy / 10))
elif [ $particle = "pi-" ]
then
    job=$((energy / 10 + 200))
fi

char="${mode:0:1}"
upper="${char^^}"
mode_upper="${upper}${mode:1}"

filename="/lustre/collider/chenjiyuan/dss-pid/run/e-signal/root/training/job${job}_${particle}_${energy}MeV/hit_sel_${particle}_${energy}MeV_ana.root"
#tree=dp
#staggered=1
event_index=235
save_dir="/lustre/collider/chenjiyuan/dss-pid/figs/"
output="${mode_upper}_${particle}_${energy}MeV.pdf"
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

#python ${mode}.py -f=$filename -i="$title" -e=$event_index -d=$save_dir -o=$output -s=$show
python ${mode}.py -f=$filename -i="$title" -e=$event_index -d=$save_dir -o="$output"