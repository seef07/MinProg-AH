# MinProg-AH

Even ter veruidelijk van waar alles ligt. Met groepje Amadeus werkte ik in de mainfolder aan een reinforcemenent Q learning algoritme. Die is nog niet af en moet nog aan geprogrammeerd worden. Die kan je vinden in /protein_folding/algorithms/reinforcement.

Ik besefte echter dat het stuctuur iets te overingewikkeld is voor zo een simpele proteine vouwprobleem als deze, dus dat verdiend ook een veel simpelere structuur. Vandaar een nieuw mapje A_new_datastru. 

## Monte-Carlo folding algoritme
te vinden in => A_new_datastru.montecarlo.py
#### algemeen
Invoer Eiwitsequentie: De sequentie van het eiwit wordt voorgesteld als een reeks karakters, waarbij elk karakter een type aminozuur vertegenwoordigt (bijvoorbeeld 'H' voor hydrofoob, 'P' voor polair, 'C' voor een ander type).

Energiematrix: Een energiematrix definieert de interactie-energie tussen paren van aminozuren. Een lagere energie wijst op gunstigere interacties.

Initiële Configuratie: Het eiwit wordt in eerste instantie uitgelegd in een rechte lijn. Elk aminozuur bezet een discrete positie op een 2D-raster.

Energieberekeningsfunctie (calculate_energy): Deze functie berekent de totale energie van een gegeven eiwitconfiguratie op basis van de energiematrix. Het neemt alleen aangrenzende aminozuren in het 2D-raster in overweging.

Validiteitscontrole (is_valid_configuration): Deze functie zorgt ervoor dat geen twee aminozuren dezelfde positie op het raster innemen.

Rotatiefunctie (rotate_segment): Deze functie draait een segment van het eiwit rond een draaipunt. De rotatie kan met de klok mee of tegen de klok in zijn. Dit is hoe het eiwit zijn configuratie verandert tijdens de simulatie.

#### Monte Carlo Vouwalgoritme (monte_carlo_folding):

Het algoritme herhaalt een gespecificeerd aantal keer.

In elke iteratie selecteert het willekeurig een draaipunt en richting om een segment van het eiwit te draaien.

Na de rotatie berekent het de energie van de nieuwe configuratie.

Als de nieuwe configuratie lagere energie heeft, wordt deze geaccepteerd.

Als de nieuwe configuratie hogere energie heeft, kan deze nog steeds worden geaccepteerd met een waarschijnlijkheid die afhangt van de temperatuur 

(die in de loop van de tijd afneemt) en het energieverschil. Dit stelt het algoritme in staat om lokale minima te vermijden.

De beste configuratie (laagste energie) gevonden tijdens de simulatie wordt bijgehouden en opgeslagen.

Visualisatiefuncties (visualize_protein, visualize_in_cmd): Deze functies helpen bij het visualiseren van de eiwitstructuur. visualize_protein gebruikt matplotlib om een grafische plot te maken, terwijl visualize_in_cmd een tekstrepresentatie print.

### WAAROM dit algoritme? 
Exploratie van Complexe Ruimten: Het Monte Carlo-algoritme is effectief in het verkennen van complexe configuratieruimten, zoals bij eiwitvouwen, waar traditionele methoden kunnen vastlopen in lokale minima.

Temperatuurafhankelijke Acceptatie: Door het accepteren van hogere-energie configuraties op basis van een afnemende 'temperatuur', kan het algoritme lokale minima vermijden en naar een globaler minimum zoeken.

### Resultaat 
Sequence : HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP

iteraties: 20000


SimulatedAnealing =  Score: -17: time: 6.3 

iterativegreed = Score: -18; time: 27.30

Monte Carlo = Score: -18; time: 9.78



Notes/TODO :
sla de bondscore elke iteratie op en plot vervolgnes de bondscore per iteratie. Vergelijk die met andere algoritmes


