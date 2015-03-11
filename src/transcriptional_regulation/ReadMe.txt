Transcriptional regulation by Pnar Pir - to reproduce Transcriptional Regulation Karr et al 2012

Python file generates the SBML file using the list of model component from the text file (tab lim)

SBML file needs two modifications before importing to Copasi: 

    <listOfUnitDefinitions>
      <unitDefinition id="substance" name="substance">
        <listOfUnits>
          <unit kind="item" exponent="1" scale="0" multiplier="1"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>

Second: Replace reversible(true) with false.

That should create intended SBML file.

Copasi file has the stochastic simulation setup and also plots the promoters, TFs and TF bound promoters. There can be only 30 promoter-TF combnations, hence 30 reactions in the model