package main

import (
	"fmt"
	"log"
	"os"

	"github.com/kpotier/selfdiff/pkg/cfg"
)

func main() {
	if len(os.Args) != 2 {
		log.Fatal("The path of the configuration file must be specified in the arguments")
	}

	log.Printf("Reading configuration file `%s`\n", os.Args[1])
	c, err := cfg.New(os.Args[1])
	if err != nil {
		log.Fatal(fmt.Errorf("newInput: %w", err))
	}

	log.Println("Calculating the velocity autocorrelation function")
	err = c.VAC()
	if err != nil {
		log.Fatal(err)
	}

	log.Println("Done")
}
