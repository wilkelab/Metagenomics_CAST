This folder contains the modified [PILER-CR](https://www.drive5.com/pilercr/) source code used for this project. 

To install, download or clone this repository and then run make:

```
git clone https://github.com/wilkelab/Metagenomics_CAST.git
cd Metagenomics_CAST/pilercr
sudo make install
```

If you do not have sudo access, you will need to manually add the pilercr executable to your $PATH. This should look something like this:

```
cd Metagenomics_CAST/pilercr
make # build the pilercr executable
mkdir ~/bin 
cp pilercr ~/bin/
```

And then add the following line to the end of your `~/.bashrc` file:

`PATH=$PATH:~/bin`
