require 'set'
require_relative 'lib/fisher_table' 
require_relative 'lib/ensembl_mappings'

raise 'Specify ranking filename'  unless ranking_filename = ARGV[0]
num_detected = ARGV[1] && ARGV[1].to_i

HgncByEnsembl = load_hgnc_by_ensembl('source_data/EnsemblToHGNC_GRCh37.p13.tsv')

protein_coding_genes = File.readlines('source_data/allhgnc_protein_coding.txt').map(&:chomp).to_set
hgncs_all = HgncByEnsembl.values.flatten.uniq.to_set # & protein_coding_genes

# oncogenes = File.readlines('source_data/сancer_genes_lists/allonco_20130923_reduced.csv').map{|l|
#   l.chomp.split
# }.map{|gene, num_evidences|
#   [gene, num_evidences.to_i]
# }

# reliable_oncogenes = oncogenes.select{|gene, num_evidences|
#   num_evidences >= 1
# }.map{|gene, num_evidences| gene}.to_set

reliable_oncogenes = File.readlines('source_data/сancer_genes_lists/onco_gathered.txt').map(&:chomp).to_set
reliable_oncogenes += ['LIPI', 'RAB11FIP3', 'TUBB8', 'SEZ6L']

known_oncogenes = reliable_oncogenes & hgncs_all

genes_ranked = File.readlines(ranking_filename)
  .map(&:chomp)
  .reject{|line| line.start_with?('#') }
  .map(&:split)
  .map{|gene, known_oncogene, num_twines, reg_length, max_mutations, uniq_mutations, max_samples, twine_maxlength, max_mutation_density, total_num_mutations, total_length|
    raise  unless ['0','1'].include?(known_oncogene)
    [gene, known_oncogene == '1' ? true : false, num_twines.to_i, reg_length.to_i, max_mutations.to_i, uniq_mutations.to_i, max_samples.to_i, twine_maxlength.to_i, max_mutation_density.to_f, total_num_mutations.to_i, total_length.to_i]
  }
  .reject{|gene, known_oncogene, num_twines, reg_length, max_mutations, uniq_mutations, max_samples, twine_maxlength, max_mutation_density, total_num_mutations, total_length|
    max_mutations <= 2 || max_samples <= 2
  }
  .sort_by{|gene, known_oncogene, num_twines, reg_length, max_mutations, uniq_mutations, max_samples, twine_maxlength, max_mutation_density, total_num_mutations, total_length|
    max_mutations
    # twine_maxlength
    # [twine_maxlength / 10.0, max_mutations].max
    # reg_length
  }
  .reverse
  .map{|gene, known_oncogene, num_twines, reg_length, max_mutations, uniq_mutations, max_samples, twine_maxlength, max_mutation_density, total_num_mutations, total_length|
    gene
  }

genes_ranked = genes_ranked.first(num_detected) if num_detected

num_detected_oncogenes = genes_ranked.inject([]){|oncogenes_hit, (gene, index)|
  oncogenes_hit_prev = oncogenes_hit.last || []
  oncogenes_hit + [oncogenes_hit_prev + (known_oncogenes.include?(gene) ? [gene] : [])]
}


puts "Known oncogenes: #{known_oncogenes.size}"
puts "All genes: #{hgncs_all.size}"

num_detected_oncogenes.each_with_index.each{|oncogenes_hit_cur, ind|
  num_predictions = ind + 1
  num_known_oncogenes_predicted = oncogenes_hit_cur.size

  ft = FisherTable.by_two_classes(
    class_a_positive: num_known_oncogenes_predicted, class_a_negative: num_predictions - num_known_oncogenes_predicted,
    class_b_positive: known_oncogenes.size - num_known_oncogenes_predicted, class_b_negative: hgncs_all.size - known_oncogenes.size - num_predictions + num_known_oncogenes_predicted,
  )

  significance = ft.significance
  # if significance < 0.05
    log_significance = significance < 0.05 ? -Math.log10(significance).round(2) : 'NO'
    puts "#{ft.a_to_b_positive_rate_ratio.round(3)}\t#{oncogenes_hit_cur.size} of #{ind + 1}\t#{log_significance}\t" + oncogenes_hit_cur.join(', ')
  # end
}
