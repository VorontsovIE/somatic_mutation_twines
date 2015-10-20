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
              }
              [chromosome, IntervalNotation::Operations.union(transcript_introns)]
            }.to_h
  result.default = IntervalNotation::Syntax::Long::Empty
  result
end

# /home/ilya/iogen/genome/hg19_exons(ensembl,GRCh37.p13).txt
def read_first_introns_by_chromosome(filename)
  result = EnsemblExon.each_in_file(filename)
            .group_by(&:chromosome)
            .map{|chromosome, exons|
              transcript_introns = exons.group_by(&:ensembl_transcript_id).map{|ensembl_transcript_id, exons|
                strand = exons.first.strand
                raise 'Different strands for the same exon list'  unless exons.all?{|exon| exon.strand == strand }
                gene_exons = IntervalNotation::Operations.union(exons.map(&:exon_region))
                all_introns = gene_exons.covering_interval - gene_exons
                # Choose only the first intron
                case strand
                when :+
                  IntervalNotation::IntervalSet.new(all_introns.intervals.first(1))
                when :-
                  IntervalNotation::IntervalSet.new(all_introns.intervals.last(1))
                else
                  raise 'Unknown strand'
                end
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


GenomeBasicInterval = Struct.new(:chromosome, :from, :to) do
  def interval
    IntervalNotation::Syntax::Long::closed_open(from + 1, to + 1) # from 0-based right excluded to 1-based right included
  end
end

def read_dhs_accessible_regions_by_chromosome(filename)
  result = File.readlines(filename).map{|line|
    chr, from, to = line.chomp.split("\t")
    GenomeBasicInterval.new(chr.sub(/^chr/, '').to_sym, from.to_i, to.to_i)
  }.group_by(&:chromosome)
  .map{|chromosome, accessible_regions|
    [chromosome, IntervalNotation::Operations.union(accessible_regions.map(&:interval))]
  }.to_h
  result.default = IntervalNotation::Syntax::Long::Empty
  result
end
