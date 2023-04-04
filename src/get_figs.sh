dir=/Volumes/TracerX/working/VHL_GERMLINE/tidda/vhl/out/230210_A01366_0351_AHNHCFDSX5_x_221202_A01366_0326_AHHTTWDMXY/hg38/comb_x_SHE5052A9_S101/all-well/DGE_filtered/unintegrated/
figs_dir=~/Desktop/figs ; mkdir $figs_dir


(
  cd $dir
  for pid in N* ; do
    echo $pid
    cp $pid/infercnv/samples/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.png \
      $figs_dir/infercnv_${pid}.png
  done
)

cp $dir/sample_heatmap_scna_arm_full_labelled.pdf $figs_dir
cp $dir/3p_loss_ridge.pdf $figs_dir
cp $dir/chr3_status.pdf $figs_dir
cp $dir/vhl_status.pdf $figs_dir
cp $dir/chr3_status_by_lesion_type.pdf $figs_dir
cp $dir/dea/normal_tme_vs_normal_normal_heatmap.pdf $figs_dir
cp $dir/vhl_status_by_partition_annot.pdf $figs_dir
cp $dir/cluster_annots_umap.pdf $figs_dir
cp $dir/sample_composition.pdf $figs_dir
