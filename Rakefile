%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require_relative "./Rakefile_common.rb"

file_base = File.expand_path(File.dirname(__FILE__)).to_s+'/lib'

task :default => [:build]

TESTS = [
  "check_1_quadratic",
  "check_2_cubic",
  "check_3_quartic"
]

desc "run tests on linux/osx"
task :run do
  TESTS.each do |cmd|
    sh "./bin/#{cmd}"
  end
end

desc "run tests (Release) on windows"
task :run_win do
  TESTS.each do |cmd|
    sh "bin\\Release\\#{cmd}.exe"
  end
end

desc "run tests (Debug) on windows"
task :run_win_debug do
  TESTS.each do |cmd|
    sh "bin\\Debug\\#{cmd}.exe"
  end
end

desc "build lib"
task :build do
  sh "make config"
  sh "make --jobs=8 install_local"
end

def ChangeOnFile( file, text_to_replace, text_to_put_in_place )
  text = File.read file
  File.open(file, 'w+'){|f| f << text.gsub(text_to_replace, text_to_put_in_place)}
end

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  FileUtils.rm_rf 'lib'

  args.with_defaults( :year => "2017", :bits => "x64" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = win_vs(args.bits,args.year)
  if COMPILE_EXECUTABLE then
    cmd_cmake += ' -DBUILD_EXECUTABLE:VAR=true '
  else
    cmd_cmake += ' -DBUILD_EXECUTABLE:VAR=false '
  end
  cmd_cmake += " -DINSTALL_HERE:VAR=true "
  #cmd_cmake += " -DCMAKE_INSTALL_PREFIX=\"#{file_base}\" "

  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/bin/"+args.bits
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"

  if COMPILE_DEBUG then
    sh cmd_cmake + ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=WARNING ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh cmd_cmake + ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=WARNING ..'
    sh 'cmake  --build . --config Release  --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc 'compile for OSX'
task :build_osx do |t, args|
  FileUtils.rm_rf 'lib'

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = 'cmake -DBUILD_EXECUTABLE:VAR='
  if COMPILE_EXECUTABLE then
    cmd_cmake += 'true '
  else
    cmd_cmake += 'false '
  end
  cmd_cmake += " -DINSTALL_HERE:VAR=true "
  #cmd_cmake += " -DCMAKE_INSTALL_PREFIX=\"#{file_base}\" "

  if COMPILE_DEBUG then
    sh cmd_cmake + '-DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=WARNING ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh cmd_cmake + ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=WARNING ..'
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end
  FileUtils.cd '..'
end

desc 'compile for LINUX'
task :build_linux do |t, args|
  FileUtils.rm_rf 'lib'

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = 'cmake -DBUILD_EXECUTABLE:VAR='
  if COMPILE_EXECUTABLE then
    cmd_cmake += 'true '
  else
    cmd_cmake += 'false '
  end
  cmd_cmake += " -DINSTALL_HERE:VAR=true "
  #cmd_cmake += " -DCMAKE_INSTALL_PREFIX=\"#{file_base}\" "

  if COMPILE_DEBUG then
    sh cmd_cmake + '-DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=WARNING ..'
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh cmd_cmake + ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=WARNING ..'
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end
  FileUtils.cd '..'
end

desc "clean for OSX"
task :clean_osx do
  FileUtils.rm_rf 'lib'
  sh "make clean"
end

desc "clean for LINUX"
task :clean_linux do
  FileUtils.rm_rf 'lib'
  sh "make clean"
end

desc "clean for WINDOWS"
task :clean_win do
  FileUtils.rm_rf 'lib'
  FileUtils.rm_rf 'vs_*'
end
