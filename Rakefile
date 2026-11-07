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

# MSYS2 non si riconosce dal solo RUBY_PLATFORM: a seconda del pacchetto ruby
# vale "x86_64-pc-msys", "x64-mingw-ucrt" o "x86_64-w64-mingw32", e in qualche
# installazione contiene "linux", nel qual caso il ramo /linux|cygwin/ scattava
# per primo e la build si dichiarava Linux pur usando il toolchain MinGW.
# MSYSTEM e' esportata da MSYS2 in ogni sua shell, quindi la si guarda prima.
$os = if !ENV['MSYSTEM'].to_s.empty? then :mingw
      else
        case RUBY_PLATFORM
        when /darwin/       then :mac
        when /msys/         then :mingw
        when /linux|cygwin/ then :linux
        else                     :win
        end
      end

$build_type = $compile_debug ? 'Debug' : 'Release'

$enable_hq = defined?(COMPILE_ARBITRARY_PRECISION) ? COMPILE_ARBITRARY_PRECISION : false

BASE_BUILD_OPTIONS = [
  "-DCMAKE_BUILD_TYPE=#{$build_type}",
  "-DQRF_BUILD_SHARED=#{$compile_dynamic ? 'ON' : 'OFF'}",
  "-DQRF_ENABLE_ARBITRARY_PRECISION=#{$enable_hq ? 'ON' : 'OFF'}"
].freeze

def cmake_options(enable_tests:)
  (BASE_BUILD_OPTIONS + [
    "-DQRF_ENABLE_TESTS=#{enable_tests ? 'ON' : 'OFF'}"
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

def configure_and_build(enable_tests: false, target: 'install')
  in_dir('build') do
    target_opt = target ? "--target #{target}" : ''
    sh "cmake -G Ninja #{cmake_options(enable_tests: enable_tests)} .."
    sh "cmake --build . --config #{$build_type} #{target_opt} #{$parallel}"
  end
end

def build_and_run_tests
  configure_and_build(enable_tests: true, target: nil)
  Dir.chdir('build') { sh "ctest -C #{$build_type} --output-on-failure" }
end

desc 'Default task: build'
task default: :build

desc 'Build with CMake/Ninja'
task :build do
  puts "Build (#{$os})".green

  configure_and_build(enable_tests: false, target: 'install')
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
  FileUtils.rm_rf(%w[build lib bin])
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
