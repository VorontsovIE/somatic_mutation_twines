require 'forwardable'
require_relative 'genome_markup/region_type'
require_relative 'sequence_with_snv'

SNVInfo = Struct.new(:variant_id, :snv_sequence,
                     :cancer_type, :sample_id, :chromosome, :position, :strand,
                     :mutation_region_types) do

  extend Forwardable
  def_delegators :mutation_region_types, *RegionType::FEATURE_INQUIRIES

  def self.from_string(str)
    variant_id, snv_sequence, \
      cancer_type, sample_id,  chromosome, position, strand, \
      mutation_region_types = str.chomp.split("\t", 8)
    new(variant_id,
        SequenceWithSNV.from_string(snv_sequence),
        cancer_type.to_sym, sample_id.to_sym,  chromosome.to_sym, position.to_i, strand.to_sym,
        RegionType.from_string(mutation_region_types))
  end

  def self.each_in_stream(stream, &block)
    stream.each_line.lazy.map{|line| from_string(line) }.each(&block)
  end

  def self.each_in_file(filename, &block)
    return enum_for(:each_in_file, filename).lazy  unless block_given?
    File.open(filename) do |f|
      f.readline # skip header
      each_in_stream(f, &block)
    end
  end

  def to_s
    SNVInfo::COLUMN_ORDER.map{|column| send(column) }.join("\t")
  end

  def marked_up(genome_markup)
    new_mutation_region_types = genome_markup.get_region_type(chromosome, position)
    SNVInfo.new(variant_id, snv_sequence,
                cancer_type, sample_id, chromosome, position, strand,
                new_mutation_region_types)
  end

  def context_before
    snv_sequence.context(before: 1, after: 1, allele_variant_number: 0)
  end

  def context_after
    snv_sequence.context(before: 1, after: 1, allele_variant_number: 1)
  end

  def reference_base
    snv_sequence.allele_variants[0]
  end

  def mutant_base
    snv_sequence.allele_variants[1]
  end

  def context_full
    snv_sequence.subsequence(before: 1, after: 1).to_s
  end

  def self.reverse_strand(original_strand)
    case original_strand
    when :+
      :-
    when :-
      :+
    else
      raise 'Unknown strand'
    end
  end

  def revcomp
    SNVInfo.new(variant_id, snv_sequence.revcomp,
                cancer_type, sample_id, chromosome, position, SNVInfo.reverse_strand(strand),
                mutation_region_types)
  end

  def in_pyrimidine_context?
    snv_sequence.in_pyrimidine_context?
  end

  def in_pyrimidine_context
    in_pyrimidine_context? ? self : revcomp
  end


  # # Deprecated
  # def load_sequence(genome_reader, five_prime_flank_length, three_prime_flank_length)
  #   genome_reader.read_sequence(chromosome, ONE_BASED_INCLUSIVE, position - five_prime_flank_length, position + three_prime_flank_length).upcase
  # end

  # # Deprecated
  # def load_site_sequence(genome_reader, site, flank_length = 0)
  #   load_sequence(genome_reader,
  #                 flank_length + site.seq_1_five_flank_length,
  #                 flank_length + site.seq_1_three_flank_length)
  # end
end

SNVInfo::COLUMN_ORDER = [:variant_id, :snv_sequence,
                         :cancer_type, :sample_id, :chromosome, :position,
                         :strand, :mutation_region_types]

SNVInfo::COLUMN_TITLES = {
  variant_id: 'Variant id', cancer_type: 'Cancer type', sample_id: 'Sample',
  chromosome: 'Chromosome', position: 'Position',
  strand: 'Strand',
  snv_sequence: 'SNV sequence',
  mutation_region_types: 'Mutation region type'
}

SNVInfo::HEADER = '# ' + SNVInfo::COLUMN_ORDER.map{|column| SNVInfo::COLUMN_TITLES[column] }.join("\t")
