# PCB Tool
Tool for modifying KiCad schematic (.kicad_sch) and PCB (.kicad_pcb) files.
Can be used to hide some or all references or remove nets from pads prior to autorouting.

## Features
* Support for KiCad 8
* Hide references using regular expression match
* Remove nets from pads using regular expression match

## Usage
| Parameter                                      | Description
|------------------------------------------------|-----------------
| -xfp \<old footprint> \<new footprint>         | Exchange footprint (only schematic)
| -xsw \<old segment width> \<new segment width> | Exchange segment width (only pcb)
| -xvs \<old via size> \<new via size>           | Exchange segment width (only pcb)
| -hr \<reference>                               | Hide reference, supports regular expressions, e.g. R.* (only pcb)
| -rn \<net>                                     | Remove net, supports regular expressions, e.g. GND (only pcb)
| -mfp \<footprint>                              | Move footprint to border (only pcb)

## Examples
| Command                                                          | Description
|------------------------------------------------------------------|-----------------
| pcb_tool -xfp R_0603_1608Metric R_0402_1005Metric test.kicad_sch | Exchange resistor footprint
| pcb_tool -hr .* test.kicad_pcb                                   | Hide all references
| pcb_tool -hr R.* test.kicad_pcb                                  | Hide all resistors
| pcb_tool -rn GND -rn \\+3V3 test.kicad_pcb                       | Remove GND and 3.3V before autorouting
