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
