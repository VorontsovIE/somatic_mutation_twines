require 'interval_notation'
require_relative 'ensembl_exon'

# /home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt
def load_promoters_by_chromosome(filename, length_5_prime: 5000, length_3_prime: 500)
  result = EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              promoter_regions = exons.select(&:transcript_start).map{|exon|
                tss_pos = exon.transcript_start
                if exon.strand == :+
                  IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_5_prime, tss_pos + length_3_prime)
                else
                  IntervalNotation::Syntax::Long.closed_closed(tss_pos - length_3_prime, tss_pos + length_5_prime)
                end
              }

              [chromosome, IntervalNotation::Operations.union(promoter_regions)]
            }.to_h
  result.default = IntervalNotation::Syntax::Long::Empty
  result
end

# /home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt
def read_introns_by_chromosome(filename)
  result = EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              transcript_introns = exons.group_by(&:ensembl_transcript_id).map{|ensembl_transcript_id, exons|
                strand = exons.first.strand
                raise 'Different strands for the same exon list'  unless exons.all?{|exon| exon.strand == strand }
                gene_exons = IntervalNotation::Operations.union(exons.map(&:exon_region))
                all_introns = gene_exons.covering_interval - gene_exons
                # # Choose only the first intron
                # case strand
                # when :+
                #   IntervalNotation::IntervalSet.new(all_introns.intervals.first(1))
                # when :-
                #   IntervalNotation::IntervalSet.new(all_introns.intervals.last(1))
                # else
                #   raise 'Unknown strand'
                # end
              }
              [chromosome, IntervalNotation::Operations.union(transcript_introns)]
            }.to_h
  result.default = IntervalNotation::Syntax::Long::Empty
  result
end

# /home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt
def read_coding_regions_by_chromosome(filename)
  result = EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              # grouping by transctripts is not necessary for now but can make sense in future
              transcript_coding_regions = exons.group_by(&:ensembl_transcript_id).map{|ensembl_transcript_id, exons|
                IntervalNotation::Operations.union(exons.map(&:coding_part_region))
              }
              [chromosome, IntervalNotation::Operations.union(transcript_coding_regions)]
            }.to_h
  result.default = IntervalNotation::Syntax::Long::Empty
  result
end



class GenomeMarkup
  attr_reader :introns_by_chromosome
  attr_reader :promoters_by_chromosome
  attr_reader :coding_regions_by_chromosome
  attr_reader :regulatory_by_chromosome

  # To load genome markup from source files, use GenomeMarkupLoader
  def initialize(introns_by_chromosome:, promoters_by_chromosome:, coding_regions_by_chromosome:)
    @introns_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(introns_by_chromosome)
    @promoters_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(promoters_by_chromosome)
    @coding_regions_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(coding_regions_by_chromosome)

    @regulatory_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty)

    chromosomes = [ @introns_by_chromosome, @promoters_by_chromosome, @coding_regions_by_chromosome ].flat_map(&:keys).uniq
    chromosomes.each do |chr|
      @regulatory_by_chromosome[chr] = ((@promoters_by_chromosome[chr] | @introns_by_chromosome[chr]) - @coding_regions_by_chromosome[chr])
    end
  end
end

# This class can store data necessary to load markup (but do it lazily and caches the result)
GenomeMarkupLoader = Struct.new(:exonic_markup_filename,
                                :promoter_length_5_prime,
                                :promoter_length_3_prime) do

  # Just a constructor with human-readable params
  def self.create(exonic_markup_filename:,
                  promoter_length_5_prime: 5000,
                  promoter_length_3_prime: 500)
    self.new(exonic_markup_filename,
            promoter_length_5_prime,
            promoter_length_3_prime)
  end

  def introns
    tm = Time.now
    @introns ||= read_introns_by_chromosome(exonic_markup_filename).tap{ $stderr.puts("introns #{Time.now - tm}") }
  end

  def promoters
    tm = Time.now
    @promoters ||= load_promoters_by_chromosome(exonic_markup_filename,
                                            length_5_prime: promoter_length_5_prime,
                                            length_3_prime: promoter_length_3_prime).tap{ $stderr.puts("promoters #{Time.now - tm}") }
  end

  def coding_regions
    tm = Time.now
    @coding_regions ||= read_coding_regions_by_chromosome(exonic_markup_filename).tap{ $stderr.puts("coding #{Time.now - tm}") }
  end

  # It can last for a very long time
  def load_markup
    @cache ||= GenomeMarkup.new(introns_by_chromosome: introns,
                               promoters_by_chromosome: promoters,
                               coding_regions_by_chromosome: coding_regions)
  end
end


