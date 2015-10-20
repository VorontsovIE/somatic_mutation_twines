require 'interval_notation'
require_relative 'load_genome_structure'
require_relative 'region_type'

class GenomeMarkup
  attr_reader :introns_by_chromosome
  attr_reader :promoters_by_chromosome
  attr_reader :coding_regions_by_chromosome
  attr_reader :regulatory_by_chromosome
  attr_reader :dhs_accessible_by_chromosome

  # To load genome markup from source files, use GenomeMarkupLoader
  def initialize(introns_by_chromosome:, promoters_by_chromosome:,
                coding_regions_by_chromosome:,
                dhs_accessible_by_chromosome: nil)
    @introns_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(introns_by_chromosome)
    @promoters_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(promoters_by_chromosome)
    @coding_regions_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(coding_regions_by_chromosome)
    if dhs_accessible_by_chromosome
      @dhs_accessible_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty).merge(dhs_accessible_by_chromosome)
    else
      @dhs_accessible_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::R)
    end

    @regulatory_by_chromosome = Hash.new(IntervalNotation::Syntax::Long::Empty)

    chromosomes = [ @introns_by_chromosome, @promoters_by_chromosome,
                    @coding_regions_by_chromosome, @dhs_accessible_by_chromosome
                  ].flat_map(&:keys).uniq
    chromosomes.each do |chr|
      @regulatory_by_chromosome[chr] = ((@promoters_by_chromosome[chr] | @introns_by_chromosome[chr]) - @coding_regions_by_chromosome[chr]) & @dhs_accessible_by_chromosome[chr]
    end
  end

  def promoter?(chromosome, position)
    case position
    when Integer
      promoters_by_chromosome[chromosome].include_position?(position)
    when IntervalNotation::IntervalSet # position is an interval
      promoters_by_chromosome[chromosome].intersect?(position)
    else
      raise TypeError, 'Position should be either integer or interval set'
    end
  end

  def intronic?(chromosome, position)
    case position
    when Integer
      introns_by_chromosome[chromosome].include_position?(position)
    when IntervalNotation::IntervalSet # position is an interval
      introns_by_chromosome[chromosome].intersect?(position)
    else
      raise TypeError, 'Position should be either integer or interval set'
    end
  end

  def coding?(chromosome, position)
    case position
    when Integer
      coding_regions_by_chromosome[chromosome].include_position?(position)
    when IntervalNotation::IntervalSet # position is an interval
      coding_regions_by_chromosome[chromosome].intersect?(position)
    else
      raise TypeError, 'Position should be either integer or interval set'
    end
  end

  def dhs_accessible?(chromosome, position)
    case position
    when Integer
      dhs_accessible_by_chromosome[chromosome].include_position?(position)
    when IntervalNotation::IntervalSet # position is an interval
      dhs_accessible_by_chromosome[chromosome].intersect?(position)
    else
      raise TypeError, 'Position should be either integer or interval set'
    end
  end

  # only singular positions
  def regulatory?(chromosome, position)
    # get_region_type(chromosome, position).regulatory?

    regulatory_by_chromosome[chromosome].include_position?(position)
  end

  def chromosome_marked_up?(chromosome)
    promoters_by_chromosome.has_key?(chromosome) || \
    introns_by_chromosome.has_key?(chromosome)
  end

  def get_region_type(chromosome, position)
    result = RegionType.new
    result << :intronic  if intronic?(chromosome, position)
    result << :promoter  if promoter?(chromosome, position)
    result << :coding    if coding?(chromosome, position)
    result << :dhs_accessible  if dhs_accessible?(chromosome, position)
    result
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
    self.new(exonic_markup_filename, promoter_length_5_prime, promoter_length_3_prime)
  end

  def introns
    @introns ||= read_introns_by_chromosome(exonic_markup_filename)
  end

  def promoters
    @promoters ||= load_promoters_by_chromosome(exonic_markup_filename,
                                            length_5_prime: promoter_length_5_prime,
                                            length_3_prime: promoter_length_3_prime)
  end

  def coding_regions
    @coding_regions ||= read_coding_regions_by_chromosome(exonic_markup_filename)
  end

  # It can last for a very long time
  def load_markup(dhs_accessible_filename: nil)
    if dhs_accessible_filename
      GenomeMarkup.new(introns_by_chromosome: introns,
                       promoters_by_chromosome: promoters,
                       coding_regions_by_chromosome: coding_regions,
                       dhs_accessible_by_chromosome: read_dhs_accessible_regions_by_chromosome(dhs_accessible_filename))
    else
      @cache ||= GenomeMarkup.new(introns_by_chromosome: introns,
                                 promoters_by_chromosome: promoters,
                                 coding_regions_by_chromosome: coding_regions,
                                 dhs_accessible_by_chromosome: nil # Hash.new(IntervalNotation::Syntax::Long::R) -- wrong way: default value will be rewritten
                                 )
    end
  end
end
