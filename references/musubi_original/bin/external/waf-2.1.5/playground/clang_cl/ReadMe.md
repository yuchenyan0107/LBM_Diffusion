# Clang-CL

Clang-CL is a drop-in MSVC compatible driver replacing CL.exe  
The clang compiler offers high compatibility with MSVC, but also offers more up to date C++ support.  
It features better code generation, but still adheres to the MSVC ABI, letting you link with link.exe, offering you superior performance but also PDB debug info.

On Windows this waf module should just work, on Linux it tries to find the LLVM replacements and requires an environment containing the paths defined by the vsvars batch files (Visual Studio C++ Developer command prompt).

# Cross compilation

To cross compile for Windows from Linux, you will require the following:

* A partition with Windows installed (NTFS).
* Visual Studio (Tested with 2017).
* The Windows SDK.
* lowntfs-3g file system driver.

Make sure the Windows partition is mounted with `-t lowntfs-3g -o defaults,ignore_case,windows_names`.  
This will allow Clang to find all headers and libraries referenced by scripts and headers, otherwise you will run into case sensitivity errors.  
You can run a script to make all filenames lowercase, but that edits your Visual Studio installation, and I don't know if that has an effect on upgradability.

Clang uses the following environment variables to detect the Visual Studio install: `VCINSTALLDIR`, `VCToolsInstallDir`, `INCLUDE`, `LIB`, `LIBPATH`  
I just copied these from the output of the `set` command in an MSVC command prompt on Windows and translated the paths to Linux paths.  
Notice how the semicolon is still used as a path separator.  
See `example_environment_linux.sh` for how my setup looks like. It expects the Windows partition to be mounted on `/mnt/windows`, with VS2017 installed and Windows 10 SDK 10.0.17763.0.

To specify a custom LLVM installation, you can put the path in the `LLVM_PATH` environment variable, or put the path in `cfg.env.LLVM_PATH` in your wscript.
