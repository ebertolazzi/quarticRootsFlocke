# frozen_string_literal: true

require 'etc'
require 'fileutils'
require 'rake/clean'

begin
  require 'colorize'
rescue LoadError
  class String
    def red
      self
    end

    def green
      self
    end

    def yellow
      self
    end
  end
end

config_files = [
  File.expand_path('../Rakefile_configure.rb', __dir__),
  File.expand_path('../../Rakefile_configure.rb', __dir__)
]

if (config = config_files.find { |path| File.exist?(path) })
  require config
end

$compile_debug = defined?(COMPILE_DEBUG) ? COMPILE_DEBUG : false
$compile_dynamic = defined?(COMPILE_DYNAMIC) ? COMPILE_DYNAMIC : false
$compile_executable = defined?(COMPILE_EXECUTABLE) ? COMPILE_EXECUTABLE : true

$os = case RUBY_PLATFORM
      when /darwin/       then :mac
      when /linux|cygwin/ then :linux
      when /msys/         then :mingw
      else                     :win
      end

$build_type = $compile_debug ? 'Debug' : 'Release'

BASE_BUILD_OPTIONS = [
  "-DCMAKE_BUILD_TYPE=#{$build_type}",
  "-DUTILS_BUILD_SHARED=#{$compile_dynamic ? 'ON' : 'OFF'}"
].freeze

def cmake_options(enable_tests:)
  (BASE_BUILD_OPTIONS + [
    "-DUTILS_ENABLE_TESTS=#{enable_tests ? 'ON' : 'OFF'}"
  ]).join(' ')
end

$parallel = if $os == :win
              ''
            else
              "--parallel #{Etc.nprocessors}"
            end

def in_dir(path)
  FileUtils.mkdir_p(path)
  Dir.chdir(path) { yield }
end

def visual_studio_arch
  cl = `where cl.exe 2>NUL`.lines.first.to_s.strip

  case cl
  when /(x64|amd64)\\cl\.exe/i then 'x64'
  when /(bin|x86|amd32)\\cl\.exe/i then 'x86'
  else
    raise 'Cannot determine Visual Studio architecture. Run from a Visual Studio Developer Prompt.'
  end
end

def configure_and_build(bits: nil, enable_tests: false, target: 'install')
  FileUtils.rm_rf('lib')
  FileUtils.rm_rf('build')

  in_dir('build') do
    bits_opt = bits ? "-DBITS=#{bits}" : ''
    target_opt = target ? "--target #{target}" : ''
    sh "cmake -G Ninja #{bits_opt} #{cmake_options(enable_tests: enable_tests)} .."
    sh "cmake --build . --config #{$build_type} #{target_opt} #{$parallel}"
  end
end

def build_and_run_tests
  bits = $os == :win ? visual_studio_arch : nil
  configure_and_build(bits: bits, enable_tests: true, target: nil)
  Dir.chdir('build') { sh "ctest -C #{$build_type} --output-on-failure" }
end

desc 'Default task: build'
task default: :build

desc 'Build with CMake/Ninja'
task :build do
  puts "Build (#{$os})".green

  bits = $os == :win ? visual_studio_arch : nil
  configure_and_build(bits: bits, enable_tests: false, target: 'install')
end

desc 'Build tests and run CTest'
task :test do
  build_and_run_tests
end

desc 'Build tests and run CTest'
task :run do
  build_and_run_tests
end

desc 'Clean build artifacts'
task :clean do
  FileUtils.rm_rf(%w[build lib])
end

desc 'Generate compile_commands.json and run cppcheck'
task :cppcheck do
  FileUtils.rm_rf('build')
  in_dir('build') do
    sh 'cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..'
    sh 'cppcheck --project=compile_commands.json'
  end
end

desc 'Run CPack from build/'
task :cpack do
  Dir.chdir('build') do
    sh 'cpack -C CPackConfig.cmake'
    sh 'cpack -C CPackSourceConfig.cmake'
  end
end

task build_osx: :build
task build_linux: :build
task build_mingw: :build
task build_win: :build
task clean_osx: :clean
task clean_linux: :clean
task clean_mingw: :clean
task clean_win: :clean
