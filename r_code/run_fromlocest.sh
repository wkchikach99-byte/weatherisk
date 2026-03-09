#! /bin/bash


if [ $# -eq 0 ]; then
    >&2 echo "No arguments provided"
    exit 1
fi

df_list=(3 5 7)
alpha_list=(07 10 13)

# Get the options
while getopts "sn:" option; do
   case $option in
      s) # display Help
         echo "Calculating for true values only";
		 df_list=(5)
		 alpha_list=(10);;
	  n) dirname=$OPTARG;;
   esac
done

if [ -z "$dirname" ]
then
  echo "Error: No directory specified (with '-n')"
  exit 1
fi
path=~/maxstable/complete/$dirname;
out=$path/logs;

# local estimates
slurmids=()
for i in "${df_list[@]}"
do
	for j in "${alpha_list[@]}"
	do
		jid2=$(sbatch -J 2_lec${i}${j} --parsable -o ${out}/2_locest_calc_${i}_${j}_%a.out c2_locest_calc.sh ${dirname} ${i} ${j})
		jid3=$(sbatch -J 3_lem${i}${j} --parsable --dependency=afterok:${jid2} -o ${out}/3_locest_merge_${i}_${j}.out c3_locest_merge.sh ${dirname} ${i} ${j})
		jid4=$(sbatch -J 4_clt${i}${j} --parsable --dependency=afterok:${jid3} -o ${out}/4_locest_clust_${i}_${j}.out c4_locest_clust.sh ${dirname} ${i} ${j})
		slurmids+=($jid4)
	done
done
# clusters

jid_clusters=$(sbatch -J 5_cltsnd --parsable --dependency=afterok:$(echo ${slurmids[*]} | tr ' ' :) -o ${out}/5_clustsnd.out c5_saunders_clust.sh ${dirname})
# inclusters
slurmids=()
for i in "${df_list[@]}"
do
	for j in "${alpha_list[@]}"
	do
		jid6=$(sbatch -J 6_inl${i}${j} --parsable --dependency=afterok:${jid_clusters} -o ${out}/6_inclusters_locest_${i}_${j}_%a.out c6_inclusters_locest.sh ${dirname} ${i} ${j})
		slurmids+=($jid6)
		jid7=$(sbatch -J 7_ilm${i}${j} --parsable --dependency=afterok:${jid6} -o ${out}/7_inclusters_locest_merge_${i}_${j}.out c7_inclusters_locest_merge.sh ${dirname} ${i} ${j})
		slurmids+=($jid7)
		jid8=$(sbatch -J 8_ins${i}${j} --parsable --dependency=afterok:${jid_clusters} -o ${out}/8_inclusters_saunders_${i}_${j}_%a.out c8_inclusters_saunders.sh ${dirname} ${i} ${j})
		slurmids+=($jid8)
		jid9=$(sbatch -J 9_ism${i}${j} --parsable --dependency=afterok:${jid8} -o ${out}/9_inclusters_saunders_merge_${i}_${j}.out c9_inclusters_saunders_merge.sh ${dirname} ${i} ${j})
		slurmids+=($jid9)
	done
done
finish=$(sbatch -J X_finish --parsable --dependency=afterok:$(echo ${slurmids[*]} | tr ' ' :) -o ${out}/X_finish.out cX_finish.sh ${dirname})