default: hyper.png

hyper.png: *.hs
	runhaskell hyper.hs

open: hyper.png
	gnome-open hyper.png
