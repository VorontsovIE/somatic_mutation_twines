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

desc 'Collect mutations and sites for Pan-cancer'
task :pan_cancer do
  pan_cancer_sites = "results/sites/PanCancer.txt"
  rm pan_cancer_sites  if File.exist?(pan_cancer_sites)
  touch pan_cancer_sites

  pan_cancer_mutations = "source_data/SNVs/PanCancer.txt"
  rm pan_cancer_mutations  if File.exist?(pan_cancer_mutations)
  touch pan_cancer_mutations

  Dir.glob('source_data/SNVs/*.txt').sort.map{|snvs_fn|
    File.basename(snvs_fn, '.txt')
  }.reject{|cancer_type|
    cancer_type == 'PanCancer'
  }.each do |cancer_type|
    `cat 'source_data/SNVs/#{cancer_type}.txt' >> #{pan_cancer_mutations}`
    `cat 'results/sites/#{cancer_type}.txt' >> #{pan_cancer_sites}`
  end
end

desc 'Collect cancer SNVs from different folder'
task :collect_SNVs do
  Dir.glob('/home/ilya/cancerSNVs_Alexandrov/results/all_introns/SNVs/Alexandrov/*').select{|fn|
    File.directory?(fn)
  }.map{|fn|
    File.basename(fn)
  }.each do |cancer_type|
    cp File.join('/home/ilya/cancerSNVs_Alexandrov/results/all_introns/SNVs/Alexandrov/', cancer_type, 'cancer.txt'), File.join('source_data/SNVs', "#{ cancer_type.gsub(/[ -]/,'') }.txt")
  end
end

directory 'source_data/motif_collection/'
directory 'source_data/motif_collection_thresholds/'
desc 'Collect motifs'
task :collect_motifs => ['source_data/motif_collection/', 'source_data/motif_collection_thresholds/']do
  `cp /home/ilya/diHOCOMOCO/final_bundle/HUMAN/mono/pwm/* source_data/motif_collection/`
  ` ls source_data/motif_collection/* | grep -P '.[S|D].pwm$' |xargs rm`
  `cp /home/ilya/diHOCOMOCO/final_bundle/HUMAN/mono/thresholds/* source_data/motif_collection_thresholds/`
  ` ls source_data/motif_collection_thresholds/* | grep -P '.[S|D].thr$' |xargs rm`
end
