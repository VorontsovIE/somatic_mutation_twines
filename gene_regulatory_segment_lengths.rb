$:.unshift File.absolute_path('./lib', __dir__)
require 'genome_markup/genome_markup'
require 'genome_markup/regulatory_markup'
require 'ensembl_mappings'

HgncByEnsembl = load_hgnc_by_ensembl('source_data/EnsemblToHGNC_GRCh37.p13.tsv')


tm = Time.now

genes_by_chromosome = load_genes_by_chromosome('source_data/ensembl_markup.tsv', length_5_prime: 5000, length_3_prime: 0)

$stderr.puts "Gene markup loaded #{Time.now - tm}"; tm = Time.now

length_by_name = genes_by_chromosome.flat_map{|chromosome, gene_by_name|
  gene_by_name.flat_map{|ensembl_name, interval|
    HgncByEnsembl[ensembl_name.to_s].map{|hgnc_name|
      [hgnc_name, interval.total_length]
    }
  }
}.group_by{|hgnc_name, length|
  hgnc_name
}.map{|hgnc_name, lengths|
  [hgnc_name, lengths.map{|_, length| length }.max]
}.reject{|hgnc_name, length| hgnc_name.nil? || hgnc_name.empty? }


File.write('gene_length.tsv', length_by_name.map{|hgnc, length| [hgnc, length].join("\t") }.join("\n"))

puts '=================='


GENOME_MARKUP_LOADER = GenomeMarkupLoader.create(
  exonic_markup_filename: 'source_data/ensembl_markup.tsv',
  promoter_length_5_prime: 5000,
  promoter_length_3_prime: 500,
)

genome_markup = GENOME_MARKUP_LOADER.load_markup

$stderr.puts "Regulatory markup loaded #{Time.now - tm}"; tm = Time.now

regulatory_length_by_name = genes_by_chromosome.flat_map{|chromosome, gene_by_name|
  regulatory = genome_markup.regulatory_by_chromosome[chromosome]
  gene_by_name.flat_map{|ensembl_name, interval|
    HgncByEnsembl[ensembl_name.to_s].map{|hgnc_name|
      [hgnc_name, (interval & regulatory).total_length]
    }
  }
}.group_by{|hgnc_name, length|
  hgnc_name
}.map{|hgnc_name, lengths|
  [hgnc_name, lengths.map{|_, length| length }.max]
}.reject{|hgnc_name, length| hgnc_name.nil? || hgnc_name.empty? }

File.write('regulatory_gene_length.tsv', regulatory_length_by_name.map{|hgnc, length| [hgnc, length].join("\t") }.join("\n"))
