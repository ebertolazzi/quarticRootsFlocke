require_relative "./cmake_utils/Rakefile_common.rb"

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  FileUtils.rm_rf 'lib'

  args.with_defaults( :year => "2017", :bits => "x64" )

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  FileUtils.mkdir_p "../lib/lib"
  FileUtils.mkdir_p "../lib/bin"
  FileUtils.mkdir_p "../lib/bin/"+args.bits
  FileUtils.mkdir_p "../lib/dll"
  FileUtils.mkdir_p "../lib/include"

  cmd_cmake = cmake_generation_command(args.bits,args.year) + cmd_cmake_build()

  puts "run CMAKE for ROOTS".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for ROOTS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc 'compile for OSX/LINUX/MINGW'
task :build_osx_linux_mingw do
  FileUtils.rm_rf 'lib'

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  cmd_cmake = "cmake " + cmd_cmake_build()

  puts "run CMAKE for ROOTS".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for ROOTS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc 'compile for LINUX'
task :build_linux do
  Rake::Task[:build_osx_linux_mingw].invoke()
end

desc 'compile for OSX'
task :build_osx do
  Rake::Task[:build_osx_linux_mingw].invoke()
end

desc 'compile for MINGW'
task :build_mingw do
  Rake::Task[:build_osx_linux_mingw].invoke()
end

desc 'pack for OSX/LINUX/MINGW/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for ROOTS".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end
