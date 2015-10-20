require_relative '../sequence_with_snv'
class SNVName
  attr_reader :name

  # PD4020a;10:100003305@A[G/A]A
  PATTERN = /\A(?<sample>.+);(?<chromosome>.+):(?<position>\d+)@(?<substitution>.+)\z/

  def initialize(name)
    @name ||= name
  end
  def to_s; name; end
  def inspect; name; end
  def matched_pattern; @match ||= PATTERN.match(name); end
  def sample_name; @sample_name ||= matched_pattern[:sample]; end
  def chromosome; @chromosome ||= matched_pattern[:chromosome].to_sym; end
  def position; @position ||= matched_pattern[:position].to_i; end
  def substitution; @substitution ||= SequenceWithSNV.from_string(matched_pattern[:substitution]); end
end
