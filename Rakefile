require_relative "./cmake_utils/Rakefile_common.rb"

task :build_common, [:bits] do |t, args|
  args.with_defaults( :bits => "" )

  puts "UTILS build (osx/linux/mingw/windows)".green

  FileUtils.rm_rf   'lib'
  FileUtils.rm_rf   'build'
  FileUtils.mkdir_p 'build'
  FileUtils.cd      'build'

  puts "run CMAKE for UTILS".yellow
  if args.bits == "" then
    sh "cmake -G Ninja " + cmd_cmake_build() + ' ..'
  else
    sh "cmake -G Ninja -DBITS:VAR=#{args.bits} " + cmd_cmake_build() + ' ..'
  end
  puts "compile with CMAKE for UTILS".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

task :mingw_pacman do
  sh 'pacman -S development'
  sh 'pacman -S mingw-w64-x86_64-toolchain'
  sh 'pacman -S mingw-w64-x86_64-cmake'
  sh 'pacman -S mingw-w64-x86_64-ninja'
end

task :build_osx   => :build_common do end
task :build_linux => :build_common do end
task :build_mingw => :build_common do end
task :build_win do
  # check architecture
  case `where cl.exe`.chop
  when /(x64|amd64)\\cl\.exe/
    VS_ARCH = 'x64'
  when /(bin|x86|amd32)\\cl\.exe/
    VS_ARCH = 'x86'
  else
    raise RuntimeError, "Cannot determine architecture for Visual Studio".red
  end
  Rake::Task[:build_common].invoke(VS_ARCH)
end

task :clean do
  FileUtils.rm_rf 'lib'
end

task :clean_osx   => :clean do end
task :clean_linux => :clean do end
task :clean_mingw => :clean do end
task :clean_win   => :clean do end

desc 'pack for OSX/LINUX/MINGW/WINDOWS'
task :cpack do
  FileUtils.cd "build"
  puts "run CPACK for ROOTS".yellow
  sh 'cpack -C CPackConfig.cmake'
  sh 'cpack -C CPackSourceConfig.cmake'
  FileUtils.cd ".."
end
