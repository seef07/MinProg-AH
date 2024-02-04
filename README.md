# MinProg-AH

<div align="center">
<figure>
    <img src="protein_folding.gif" height="420">
    <h4></h4>
    <figcaption>A sample 3D protein fold.</figcaption>
</figure>
</div>

## How to Use

The Protein Folding Simulator allows you to simulate protein folding using various algorithms. Follow these instructions to get started:

1. **Installation**

   Make sure you have Python installed on your system.

2. **Clone the Repository**

   Clone this repository to your local machine

   Make sure to install the requirements that oyu can find in requirements.txt


You will be prompted to enter an algorithm choice (an integer from 1 to 7) for simulating protein folding.

5. **Choose an Algorithm**

Choose one of the available algorithms by entering the corresponding integer:

- 1: Monte Carlo
- 2: Monte Carlo 3D
- 3: Bruteforce
- 4: Greedy
- 5: Iterative Greedy
- 6: Iterative Random
- 7: Simulated Annealing

6. **Follow the On-Screen Instructions**

Follow the on-screen instructions to input your sequence and algorithm choice.

7. **View Results**

Depending on your algorithm choice, you will see different results and visualizations related to protein folding. The output is for every algorthm different. Make sure to take a look at the code to underrstand what you see. 

8. **Experiment and Explore**

Feel free to explore and experiment with different algorithms.

Have fun simulating protein folding with the Protein Folding Simulator!

## GENERAL
Protein Folding Environment
Protein Sequence: The protein sequence is the linear sequence of amino acids represented by characters such as 'H' (Hydrophobic), 'P' (Polar), and 'C' (Cysteine). This sequence is provided as input to the algorithm.

Energy Matrix: An energy matrix defines the interaction energies between pairs of amino acids. Different pairs may have attractive or repulsive interactions, affecting the overall energy of the protein configuration.

## Algoritms

### Reinforcement Learning Approach for Protein Folding
This is a completely self developed learning-based approach for protein folding. This approach aims to learn a policy that guides the folding process by selecting actions (rotations) for each amino acid in the protein sequence. 

#### Reinforcement Learning Components
##### Policy Weights: The algorithm maintains a set of policy weights, represented as a list. These weights determine the preference for different actions (rotations) based on the estimated rewards.

##### Action Selection: 
The algorithm uses the policy weights to select actions for each amino acid. The action selection process takes into account both the policy weights and a temperature parameter that controls the exploration-exploitation trade-off.

##### Episode:
 An episode in reinforcement learning corresponds to one folding attempt. During each episode, the algorithm aims to improve the folding of the protein.

##### Reward Function:
 The reward function computes the energy of the protein configuration. Lower energy values represent more stable protein conformations, and the goal is to minimize energy.

##### Heuristics: 
Several heuristics are used to estimate the quality of different actions. These heuristics include compactness, folding patterns, and distance measures. The heuristics contribute to the total estimated reward and influence action selection.

##### Learning Rate:
 The learning rate parameter controls how much the policy weights are updated based on the rewards received during the episode.

### Reinforcement Learning Process
#### Initialization:
 The algorithm initializes the state, which represents the current positions of amino acids in the protein sequence, and sets up the policy weights.

#### Episode Execution: 
The algorithm runs multiple episodes (folding attempts). Each episode consists of the following steps:

#### Action Selection:
 For each amino acid, the algorithm selects an action (rotation) based on the current state, policy weights, and temperature. The selected actions may involve clockwise or counterclockwise rotations.

#### State Transition:
 The selected actions are applied to the current state, resulting in a new state. The algorithm computes the reward for this new state based on the energy of the protein configuration.

#### Policy Weight Update:
 The policy weights are updated based on the rewards received during the episode. The update considers the influence of heuristics and the learning rate.

#### Exploration-Exploitation Trade-off: !!!! 
The algorithm balances exploration (trying new actions) and exploitation (selecting actions based on policy weights) using the temperature parameter. As episodes progress, the temperature decreases, making the algorithm more selective in action choices.

##### Randomizing Weights: 
Periodically, the algorithm randomizes the policy weights to introduce exploration and prevent convergence to suboptimal solutions.

### MonteCarlo Algoritm

#### Energy Matrix: 
Each pair of amino acids interacts differently. An "energy matrix" defines the energy associated with interactions between different amino acid pairs. For example, some pairs may attract each other (negative energy), while others may repel (positive energy).

#### Energy Calculation: 
To assess a protein's conformation (3D arrangement), we calculate its total energy based on the positions of amino acids and the sequence. The energy considers interactions between pairs of amino acids and their spatial proximity.

#### Configuration Validity Check: 
We ensure that the protein configuration is valid, meaning there are no overlaps or steric clashes between amino acids. Overlapping amino acids would represent an invalid, physically unrealistic structure.

#### Segment Rotation: 
We simulate the folding process by rotating segments of the protein. These segments consist of consecutive amino acids. We can rotate a segment clockwise or counterclockwise to explore different conformations.

#### Monte Carlo Folding:
 We use the Monte Carlo method, which involves repeated random sampling, to simulate the folding process. The simulation is based on a temperature parameter and energy differences. At higher temperatures, the algorithm can accept conformational changes even if they increase energy. As the temperature decreases, the algorithm becomes more selective in accepting higher-energy conformations, leading to stable structures.

#### Visualization:
 To understand the folding progress, we can visualize the protein configuration at different stages. This helps us see how the protein evolves from an extended chain to a compact, stable structure.

#### Bond Score History:
 We track the bond scores (energy) throughout the simulation iterations. This history allows us to observe how the algorithm improves the protein's conformation over time.

## General Algorithms exlainned

### BruteForce

This algorithm attempts to test out every imaginable way of folding a protein. In theory, given an unlimited amount of time, it could eventually find the absolute best configuration. However, in practice, it becomes highly impractical and time-consuming, especially for larger and more complex proteins.

### Iterative Random

Iterates through the protein, examining the potential movement directions for each individual node while taking precautions to avoid any collisions.

### Greedy 

Traverses the protein structure, flexing each node towards the most advantageous orientation determined by its score or, alternatively, selects a random direction when no clear preference exists.

### Greedy Iterative

Iterates through the protein, optimizing node positions by consistently selecting the most advantageous direction based on the score, steadily improving the overall protein configuration.

### Regression

Makes a number of random folds, only keeping the new shape if it has a better score than the previous shape.

### Simulated Annealing

 Resembles the Regression method, with the distinction that there is a gradually diminishing probability, in each iteration, of accepting a new protein shape, even if it happens to be less favorable than the current one.

### Spiral

Arranges the protein in a spiral shape

###
#### Heuristieken:
Het algoritme maakt gebruik van vijf heuristieken om beslissingen te nemen over hoe de aminozuren moeten worden verplaatst. Deze heuristieken omvatten zaken als het proberen te voorkomen van ongunstige interacties tussen aminozuren, het maximaliseren van compactheid en het volgen van een patroon in de sequentie.

#### Vouwen:
In elke iteratie van het algoritme wordt een aminozuur geselecteerd op basis van de heuristieken, en vervolgens wordt bepaald of het aminozuur met de klok mee of tegen de klok in moet worden verplaatst. Dit wordt gedaan om de configuratie te optimaliseren volgens de heuristieken en de energiematrix.

#### Beloning en Evaluatie:
Nadat een aminozuur is verplaatst, wordt de resulterende configuratie geëvalueerd aan de hand van een beloningsfunctie. Deze functie meet hoe goed de nieuwe configuratie is in termen van energie, waarbij gunstige configuraties lagere energiewaarden hebben.

#### Update:
De algoritme past vervolgens zijn beleid aan op basis van de beloningen die zijn ontvangen. Als een actie positief heeft bijgedragen aan het verminderen van de energie, worden de gewichten van de bijbehorende heuristieken verhoogd. Als een actie negatief heeft bijgedragen, worden de gewichten verlaagd.

#### Herhaling:
Het bovenstaande proces wordt herhaald voor meerdere iteraties of "episodes". Gedurende deze afleveringen wordt geprobeerd om de energie van de configuratie te minimaliseren en zo een optimale vouwing te bereiken.

### Waarom dit algoritme? 
Eenvoudig aan te passen: Het algoritme maakt gebruik van heuristieken en gewichten, waardoor het gemakkelijk kan worden aangepast aan verschillende eiwitten of problemen door de heuristieken en gewichten aan te passen.
Lokale zoekmethode: Het is een lokale zoekmethode, wat betekent dat het probeert de huidige oplossing te verbeteren zonder het gehele zoekruimte te verkennen, wat efficiënt kan zijn voor grote eiwitsequenties.
Flexibiliteit: Het algoritme kan worden aangepast om verschillende doelstellingen en beperkingen in te bouwen, afhankelijk van de specifieke behoeften van het probleem.
Heuristische aanpak: Het maakt gebruik van heuristieken die zijn gebaseerd op intuïtie en kennis van het eiwitten vouwprobleem, wat kan helpen bij het vinden van goede oplossingen.

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

Run1 : Reinforcement learning: Score: -21  [[0, 0, 'H'], [1, 0, 'H'], [2, 0, 'P'], [2, 1, 'C'], [1, 1, 'H'], [1, 2, 'H'], [1, 3, 'P'], [2, 3, 'C'], [2, 2, 'C'], [3, 2, 'P'], [3, 3, 'C'], [4, 3, 'P'], [5, 3, 'P'], [6, 3, 'H'], [6, 4, 'H'], [6, 5, 'H'], [5, 5, 'H'], [5, 4, 'P'], [4, 4, 'P'], [3, 4, 'H'], [2, 4, 'C'], [1, 4, 'H'], [0, 4, 'P'], [-1, 4, 'H'], [-1, 3, 'P'], [0, 3, 'H'], [0, 2, 'C'], [0, 1, 'H'], [-1, 1, 'P'], [-1, 0, 'P']]

Run 2: Reinforcement learning score: -22 ->: [[0, 0, 'H'], [1, 0, 'H'], [2, 0, 'P'], [2, -1, 'C'], [1, -1, 'H'], [0, -1, 'H'], [0, -2, 'P'], [1, -2, 'C'], [2, -2, 'C'], [2, -3, 'P'], [1, -3, 'C'], [1, -4, 'P'], [1, -5, 'P'], [2, -5, 'H'], [3, -5, 'H'], [4, -5, 'H'], [4, -4, 'H'], [3, 
-4, 'P'], [3, -3, 'P'], [3, -2, 'H'], [3, -1, 'C'], [4, -1, 'H'], [4, -2, 'P'], [4, -3, 'H'], [5, -3, 'P'], [5, -4, 'H'], [5, -5, 'C'], [6, -5, 'H'], [7, -5, 'P'], [8, -5, 'P']]
weights run 2: [-0.22899999999999998, -0.10699999999999998, 8.84000000000054, -0.07, 10.988999999999349]

Run 3: Reinforcement leanring: bestscore: -24 Beststate: [[0, 0, 'H'], [1, 0, 'H'], [1, 1, 'P'], [2, 1, 'C'], [3, 1, 'H'], [4, 1, 'H'], [4, 0, 'P'], [3, 0, 'C'], [2, 0, 'C'], [2, -1, 'P'], [3, -1, 'C'], [3, -2, 'P'], [2, -2, 'P'], [1, -2, 'H'], [1, -1, 'H'], [0, -1, 'H'], [0, -2, 'H'], [0, -3, 'P'], [-1, -3, 'P'], [-1, -2, 'H'], [-1, -1, 'C'], [-1, 0, 'H'], [-1, 1, 'P'], [-1, 2, 'H'], [0, 2, 'P'], [1, 2, 'H'], [2, 2, 'C'], [3, 2, 'H'], [3, 3, 'P'], [4, 3, 'P']]
[-0.26, -0.12, 1.088999999999991, -0.07, 2.0199999999998886]


#### Langere sequence ?
Sequence: PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP

Iterative greedy. Score: -5 ; time: 30.081570625305176 seconds

simulatedAnealing. Score: -9; time: 6.57  

Monte carlo. Score: -11; Time: Unknown but not more then 10 seconds.


Notes/TODO :
sla de bondscore elke iteratie op en plot vervolgnes de bondscore per iteratie. Vergelijk die met andere algoritmes
Convergence met temperatuur visualiseren
Opdrachten: Objective function -> UpperBound Lowerbound; Constraints toevoegen. Minimaal 2.
Introductie - methode - resultaten - discussie - conclusie
Fress
