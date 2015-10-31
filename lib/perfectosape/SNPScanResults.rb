require 'set'
require 'forwardable'
require 'interval_notation'

module PerfectosAPE
  Result = Struct.new(:variant_id, :motif_name,
                      :fold_change, :pvalue_1, :pvalue_2,
                      :pos_1, :pos_2, :length
                      ) do

    # "27610826_3 MAZ_f1  -7  direct  cggctgaGgaggaggag -7  direct  cggctgaCgaggaggag G/C 1.1218764110455249E-4 9.602413003842941E-4  0.11683275970285215"
    def self.from_string(line)
      variant_id, motif_name,
                pos_1, _orientation_1, seq_1,
                pos_2, _orientation_2, _seq_2,
                _variants,
                pvalue_1, pvalue_2, fold_change = line.split("\t")
      self.new( variant_id, motif_name.to_sym,
                fold_change.to_f, pvalue_1.to_f, pvalue_2.to_f,
                pos_1.to_i, pos_2.to_i, seq_1.length
              )
    end

    def to_s
      [variant_id, motif_name,
      pos_1, pos_2, length,
      pvalue_1, pvalue_2, fold_change].join("\t")
    end

    def site_before_substitution?(pvalue_cutoff: 0.0005)
      pvalue_1 <= pvalue_cutoff
    end

    def site_after_substitution?(pvalue_cutoff: 0.0005)
      pvalue_2 <= pvalue_cutoff
    end

    def disrupted?(fold_change_cutoff: 5)
      fold_change <= (1.0 / fold_change_cutoff)
    end

    def emerged?(fold_change_cutoff: 5)
      fold_change >= fold_change_cutoff
    end

    def self.each_in_stream(stream, &block)
      stream.each_line.lazy.reject{|line|
        line.start_with?('#')
      }.map{|line|
        self.from_string(line)
      }.each(&block)
    end

    def self.each_in_file(all_mutations_filename, &block)
      return enum_for(:each_in_file, all_mutations_filename).lazy  unless block_given?
      File.open(all_mutations_filename) do |f|
        each_in_stream(f, &block)
      end
    end


    # Task-specific methods
    require_relative 'snv_name'
    def snv_name
      @snv_name ||= SNVName.new(variant_id)
    end

    extend Forwardable
    def_delegators :snv_name, :chromosome, :sample_name, :substitution
    def_delegator :snv_name, :position, :mutation_genomic_position

    # Always takes the strongest site among two alleles
    def strongest_site_interval
      if fold_change <= 1 # affinity loss
        IntervalNotation::Syntax::Long.closed_closed(mutation_genomic_position + pos_1, mutation_genomic_position + pos_1 + length)
      else  # affinity gain
        IntervalNotation::Syntax::Long.closed_closed(mutation_genomic_position + pos_2, mutation_genomic_position + pos_2 + length)
      end
    end
  end
end
