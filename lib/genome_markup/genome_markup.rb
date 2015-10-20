require 'interval_notation'
require_relative 'ensembl_exon'

# closed boundaries; 1-based notation
def region_by_boundaries(left, right)
  if left && right
    if left != right
      IntervalNotation::Syntax::Long.closed_closed(left, right)
    else
      IntervalNotation::Syntax::Long.point(left)
    end
  else
    if !left && !right
      IntervalNotation::Syntax::Long::Empty
    else
      raise 'Wrong parsing assumptions. My bug!'
    end
  end
end

def take_one_of_same(collection, &block)
  results = collection.map(&block)
  uniq_results = results.uniq
  raise "Elements should be the same: #{results.join('; ')}"  if uniq_results.size > 1
  uniq_results.first
end

def region_around_position(pos, strand, length_5_prime:, length_3_prime:)
  if strand == :+
    region_by_boundaries(pos - length_5_prime, pos + length_3_prime)
  elsif strand == :-
    region_by_boundaries(pos - length_3_prime, pos + length_5_prime)
  else
    raise 'Unknown strand'
  end
end

def region_around_region(pos_from, pos_to, strand, length_5_prime:, length_3_prime:)
  if strand == :+
    region_by_boundaries(pos_from - length_5_prime, pos_to + length_3_prime)
  elsif strand == :-
    region_by_boundaries(pos_from - length_3_prime, pos_to + length_5_prime)
  else
    raise 'Unknown strand'
  end
end


def transcript_promoter(transcript_exons, length_5_prime: 5000, length_3_prime: 500)
  strand = take_one_of_same(transcript_exons, &:strand)
  tss_pos = take_one_of_same(transcript_exons, &:transcript_start)
  region_around_position(tss_pos, strand, length_5_prime: length_5_prime, length_3_prime: length_3_prime)
  promoter_region
end

def load_promoters_by_chromosome(filename, length_5_prime: 5000, length_3_prime: 500)
  EnsemblExon.each_in_file(filename)
  .group_by(&:chromosome)
  .map{|chromosome, exons|
    promoters_by_gene = exons.group_by(&:ensembl_gene_id).map{|ensembl_gene_id, gene_exons|
      gene_promoters = gene_exons
        .group_by(&:ensembl_transcript_id)
        .map{|ensembl_transcript_id, transript_exons|
          transcript_promoter(transript_exons, length_5_prime: length_5_prime, length_3_prime: length_3_prime)
        }
      
      promoters_region = IntervalNotation::Operations.union(gene_promoters)
      [ensembl_gene_id, promoters_region]
    }.to_h
    [chromosome, promoters_by_gene]
  }.to_h
end

def load_genes_by_chromosome(filename, length_5_prime: 5000, length_3_prime: 0)
  EnsemblExon.each_in_file(filename)
  .group_by(&:chromosome)
  .map{|chromosome, exons|
    gene_regions_by_name = exons.group_by(&:ensembl_gene_id).map{|ensembl_gene_id, gene_exons|
      gene_start = take_one_of_same(gene_exons, &:gene_start)
      gene_end = take_one_of_same(gene_exons, &:gene_end)
      strand = take_one_of_same(gene_exons, &:strand)
      whole_gene_region = region_around_region(gene_start, gene_end, strand, length_5_prime: length_5_prime, length_3_prime: length_3_prime)
      [ensembl_gene_id, whole_gene_region]
    }.to_h
    [chromosome, gene_regions_by_name]
  }.to_h
end


def gene_introns(gene_exons, length_5_prime: 5000, length_3_prime: 500)
  gene_start = take_one_of_same(gene_exons, &:gene_start)
  gene_end = take_one_of_same(gene_exons, &:gene_end)
  whole_gene_region = region_by_boundaries(gene_start, gene_end)
  exonic_region = IntervalNotation::Operations.union(gene_exons.map(&:exon_region))
  whole_gene_region - exonic_region
end

def load_introns_by_chromosome(filename, length_5_prime: 5000, length_3_prime: 500)
  EnsemblExon.each_in_file(filename)
  .group_by(&:chromosome)
  .map{|chromosome, exons|
    introns_by_gene = exons.group_by(&:ensembl_gene_id).map{|ensembl_gene_id, gene_exons|
      introns = gene_introns(gene_exons, length_5_prime: length_5_prime, length_3_prime: length_3_prime)
      [ensembl_gene_id, introns]
    }.to_h
    [chromosome, introns_by_gene]
  }.to_h
end
