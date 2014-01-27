require "preflight"

preflight = Preflight::Profiles::PDFA1A.new

puts preflight.check("thesis.pdf").inspect

File.open("thesis.pdf", "rb") do |file|
  puts preflight.check(file).inspect
end
