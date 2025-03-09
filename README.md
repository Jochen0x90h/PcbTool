# PCB Tool

Tool for modifying KiCad PCB (.kicad_pcb) files.
Can be used to hide some or all references or remove nets from pats prior to autorouting.

Examples:

| Description                            | Command
|----------------------------------------|-----------------
| Hide all references                    | pcb_tool -hr .*
| Remove GND and 3.3V before autorouting | pcb_tool -rn GND -rn \\+3V3

## Features
* Support for KiCad 8
* Hide references using regular expression match
* Remove nets from pads using regular expression match

## Usage
Pass a .kicad_pcb file as first argument. The tool generates BOM.csv and CPL.csv files next to it
