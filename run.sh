GH_USER="metapozza"; AUTHOR_NAME="Matteo Pozza"; AUTHOR_EMAIL="metapozza@proton.me"

sed -i.bak "s|^readme = \"README.md\"|readme = \"README.md\"\nauthors = [\"$AUTHOR_NAME <$AUTHOR_EMAIL>\"]\nrepository = \"https://github.com/$GH_USER/sigring\"\nhomepage = \"https://github.com/$GH_USER/sigring\"\ndocumentation = \"https://docs.rs/sigring\"|" Cargo.toml
sed -i.bak "s|AUTHOR_NAME|$AUTHOR_NAME|; s|GH_USER|$GH_USER|" NOTICE
sed -i.bak "s|AUTHOR_EMAIL|$AUTHOR_EMAIL|" COMMERCIAL.md
rm -f *.bak
