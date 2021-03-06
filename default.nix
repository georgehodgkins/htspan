with import <nixpkgs> {};
stdenv.mkDerivation rec {
	name = "env";
	env = buildEnv { name = name; paths = buildInputs; };
	buildInputs = [ gsl lzma bzip2 atlas ];
}
