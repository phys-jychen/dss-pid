#! /bin/bash
source /home/chenjiyuan/conda.env

# Available modes: ‘display’ for event display, and ‘projection’ for energy projection.  They must match the names of the existing python files!
# Available event types:  signal;  en_ecal, en_target, gmm_ecal, gmm_target, pn_target;  inclusive.
mode="display"
event_type="inclusive"
mass=5
job=5
run_index=0
event_index=2000999

char="${mode:0:1}"
upper="${char^^}"
mode_upper="${upper}${mode:1}"

filename="/lustre/collider/chenjiyuan/dss-pid/run/dp-signal/"
output="${mode_upper}_${event_type}"
if [ ${event_type} = "signal" ]
then
    filename+="${event_type}/root/test/Mass${mass}MeV/hit_Mass${mass}MeV.root"
    output+="_${mass}MeV.pdf"
elif [ ${event_type} = "inclusive" ]
then
    filename+="${event_type}/root/test/job${job}/hit_${event_type}_${job}.root"
    output+=".pdf"
elif [ ${event_type} = "en_ecal" ] || [ ${event_type} = "en_target" ] || [ ${event_type} = "gmm_ecal" ] || [ ${event_type} = "gmm_target" ] || [ ${event_type} = "pn_target" ]
then
    filename+="rare/${event_type}/test/job${job}/hit_${event_type}_${job}.root"
    output+=".pdf"
fi
#tree=dp
#staggered=1
save_dir="/lustre/collider/chenjiyuan/dss-pid/figs/"
#show=0

title=""
if [ ${event_type} = "signal" ]
then
    title+="${mass} MeV Signal"
elif [ ${event_type} = "inclusive" ]
then
    title+="Inclusive Background"
elif [ ${event_type} = "en_ecal" ]
then
    title+="EN ECAL"
elif [ ${event_type} = "en_target" ]
then
    title+="EN Target"
elif [ ${event_type} = "gmm_ecal" ]
then
    title+="GMM ECAL"
elif [ ${event_type} = "gmm_target" ]
then
    title+="GMM Target"
elif [ ${event_type} = "pn_target" ]
then
    title+="PN Target"
fi

#python ${mode}.py -f=${filename} -i="${title}" -r=${run_index} -e=${event_index} -d=${save_dir} -o=${output} -s=${show}
python ${mode}.py -f=${filename} -i="${title}" -r=${run_index} -e=${event_index} -d=${save_dir} -o=${output}