require 'perfectosape/SNPScanRunner'
desc 'Collect sites lying over SNVs'
task :collect_sites

Dir.glob('source_data/SNVs/*.txt').sort.each do |snvs_fn|
  cancer_type = File.basename(snvs_fn, File.extname(snvs_fn))

  task "collect_sites:#{cancer_type}" do
    mkdir_p 'results/sites/'
    output_file = "results/sites/#{cancer_type}.txt"

    Ape.run_SNPScan(
      snvs_file: snvs_fn,
      motif_collection: 'source_data/motif_collection',
      precalulated_thresholds: 'source_data/motif_collection_thresholds',
      fold_change_cutoff: Configuration::FoldChangeCutoff,
      pvalue_cutoff: Configuration::PvalueCutoff,
      jvm_options: ['-Xmx1G'],
      output_file: output_file
    )
  end
  task :collect_sites => "collect_sites:#{cancer_type}"
end
