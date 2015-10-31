directory 'results/twines'

desc 'Make twines'
task 'aggregate_twines' => 'results/twines'

Dir.glob('results/sites/*.txt').map{|fn|
  File.basename(fn, '.txt')
}.each do |cancer_type|
  task "aggregate_twines:#{cancer_type}" => 'results/twines/' do
    sh 'ruby', 'glue_intervals.rb', cancer_type, out: "results/twines/#{cancer_type}.yaml"
  end
  task 'aggregate_twines' => "aggregate_twines:#{cancer_type}"
end
