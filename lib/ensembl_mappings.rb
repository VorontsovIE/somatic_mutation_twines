EnsemblMappings = Struct.new(:ensg, :enst, :hgnc_symbol, :description, :uniprot_id, :uniprot_ac) do
  def self.from_string(str)
    ensg, enst, hgnc_symbol, description, uniprot_id, uniprot_ac = str.chomp.split("\t", 6)
    self.new(ensg, enst, hgnc_symbol, description, uniprot_id, uniprot_ac)
  end

  def self.each_in_file(filename, &block)
    File.readlines(filename).map{|line|
      self.from_string(line)
    }.each(&block)
  end
end

def load_hgnc_by_ensembl(mapping_filename)
  result = EnsemblMappings.each_in_file(mapping_filename).group_by(&:ensg).map{|ensg, infos|
    [ensg, infos.map(&:hgnc_symbol).uniq.compact]
  }.to_h
  result.default_proc = ->(h,k){ h[k] = [] }
  result
end
