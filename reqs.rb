$:.unshift File.absolute_path('./lib', __dir__)
require_relative 'configuration'
require 'set'
require 'snv_info'
require 'perfectosape/SNPScanResults'
require 'interval_notation'
require 'yaml'
require 'ensembl_mappings'
require 'twine'
