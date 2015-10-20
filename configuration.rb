$:.unshift File.absolute_path('./lib', __dir__)
require 'genome_markup/genome_markup'
module Configuration
  PvalueCutoff = 0.0005
  FoldChangeCutoff = 4
end
