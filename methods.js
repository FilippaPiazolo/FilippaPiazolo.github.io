/*
Settings stored in localStorage:
Genom Type: #genomeType
Genome Reference: #genomeRef
Data Type 1: #datatype1
Data Type 2: #datatype2

Dataset 1: #data1 
Dataset 2: #data2
Datafile1: #filename1
Datafile2: #filename2
Number of mutations dataset 1: #numdata1
Number of mutations dataset 2: #numdata2
Chromosome Number: #chrNumber

Storing Parameters:
only Density:
Bin Size: #binSize1
Bandwith: #bandwith1
Amplitude: #amplitude1
Min. Number of Mutations: #numMut
Thresholdvalue for cluster: #cluster
Min. Cluster Size: #minClusterSize

Comparison Density: 
Bin Size: #binSize2
Bandwith: #bandwith2
Amplitude: #amplitude2
*/


// Function for generating a text file URL containing given text
function generateTextFileUrl(txt) {
    let fileData = new Blob([txt], {type: 'text/plain'});

    var textFileUrl = window.URL.createObjectURL(fileData);

    // Returns a reference to the global variable holding the URL
    // Again, this is better than generating and returning the URL itself from the function as it will eat memory if the file contents are large or regularly changing
    return textFileUrl;
};

function downloadSVG(svgElement,downloadLink){
  var svg = document.getElementById(svgElement);
  svg.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  svg.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  svg = document.getElementById(svgElement).outerHTML;
  var svgBlob = new Blob([svg], {type: "image/svg"});
  document.getElementById(downloadLink).href = window.URL.createObjectURL(svgBlob);
}

//--------------- defined values ------------------ //
const human_genome = [
  {ID: '1', size: '249239465', centromere: '121535434,124535434'}, 
  {ID: '2', size: '243199373', centromere: '92326171,95326171'},
  {ID: '3', size: '199411731', centromere: '90504854,93504854'},
  {ID: '4', size: '191252270', centromere: '49660117,52660117'},
  {ID: '5', size: '180915260', centromere: '46405641,49405641'},
  {ID: '6', size: '171115067', centromere: '58830166,61830166'},
  {ID: '7', size: '159138663', centromere: '58054331,61054331'},
  {ID: '8', size: '146364022', centromere: '43838887,46838887'},
  {ID: '9', size: '141213431', centromere: '47367679,50367679'},
  {ID: '10', size: '135534747', centromere: '39254935,42254935'},
  {ID: '11', size: '135006516', centromere: '51644205,54644205'},
  {ID: '12', size: '133851895', centromere: '34856694,37856694'},
  {ID: '13', size: '115169878', centromere: '16000000,19000000'},
  {ID: '14', size: '107349540', centromere: '16000000,19000000'},
  {ID: '15', size: '102531392', centromere: '17000000,20000000'},
  {ID: '16', size: '90354753', centromere: '35335801,38335801'},
  {ID: '17', size: '81195210', centromere: '22263006,25263006'},
  {ID: '18', size: '78077248', centromere: '15460898,18460898'},
  {ID: '19', size: '64705560', centromere: '24681782,27681782'},
  {ID: '20', size: '63025520', centromere: '26369569,29369569'},
  {ID: '21', size: '48129895', centromere: '11288129,14288129'},
  {ID: '22', size: '51304566', centromere: '13000000,16000000'},
  {ID: 'X', size: '155270560', centromere: '58632012,61632012'},
  {ID: 'Y', size: '59373566', centromere: '10104553,13104553'}
];

const mouse_genome = [
  {ID: "1", size: "195471971", centromere: "110000,3000000"},
  {ID: "2", size: "182113224", centromere: "110000,3000000"},
  {ID: "3", size: "160039680", centromere: "110000,3000000"},
  {ID: "4", size: "156508116", centromere: "110000,3000000"},
  {ID: "5", size: "151834684", centromere: "110000,3000000"},
  {ID: "6", size: "149736546", centromere: "110000,3000000"},
  {ID: "7", size: "145441459", centromere: "110000,3000000"},
  {ID: "8", size: "129401213", centromere: "110000,3000000"},
  {ID: "9", size: "124595110", centromere: "110000,3000000"},
  {ID: "10", size: "130694993", centromere: "110000,3000000"},
  {ID: "11", size: "122082543", centromere: "110000,3000000"},
  {ID: "12", size: "120129022", centromere: "110000,3000000"},
  {ID: "13", size: "120421639", centromere: "110000,3000000"},
  {ID: "14", size: "124902244", centromere: "110000,3000000"},
  {ID: "15", size: "104043685", centromere: "110000,3000000"},
  {ID: "16", size: "98207768", centromere: "110000,3000000"},
  {ID: "17", size: "94987271", centromere: "110000,3000000"},
  {ID: "18", size: "90702639", centromere: "110000,3000000"},
  {ID: "19", size: "61431566", centromere: "110000,3000000"},
  {ID: "X", size: "171031299", centromere: "110000,3000000"},
  {ID: "Y", size: "91744698", centromere: "110000,3000000"}
]

var human_centromere = human_genome.map(human_genome => human_genome.centromere.split(/,/));
var mouse_centromere = mouse_genome.map(mouse_genome => mouse_genome.centromere.split(/,/));


var color = d3.schemePaired;
var color1 = ['DarkOrange', 'green', 'red', 'purple', 'pink', 'blue']//['#1F77B4FF', '#FF7F0EFF', '#2CA02CFF', '#D62728FF', '#9467BDFF', '#E377C2FF'];//color.slice(0,6);//["green", "blue", "red", "purple", "pink", "orange"];
var color2 = ['turquoise', 'LawnGreen', 'gold', 'SkyBlue', 'Coral', 'BurlyWood']//['#17BECFFF', '#B5CF6BFF', '#E7BA52FF', '#FF9896FF', '#C49C94FF', '#9C9EDEFF'];//color.slice(6);



//---------------- CONVERT DATASET INTO RIGHT FORMAT -------------------//
//convert chromosome number to only the number (from chr1 --> 1; chrX --> 23/20)
function convert_data(data, split) { 
  var dataset = [];
  var samples = [];

  var text = data.replaceAll("\"","");
  var ar = text.split("\n");

  var chr = localStorage.getItem("h_chr");
  var pos = localStorage.getItem("h_pos");
  var start = localStorage.getItem("h_start");
  var end = localStorage.getItem("h_end");
  var gene = localStorage.getItem("h_gene");
  var sample = localStorage.getItem("h_sample");
  var value = localStorage.getItem("h_value"); 

  //convert seperator
  if (split == 'tab')
    split = '\t';
  else if (split == "comma")
    split = ',';
  else if (split == "semicolon")
    split = ';'; 
 

  var def = ar[0].split(split);
  var sampleIndex = def.indexOf(sample);
  var chrIndex = def.indexOf(chr);
  var posIndex = def.indexOf(pos);
  var startIndex = def.indexOf(start);
  var endIndex = def.indexOf(end);
  var geneIndex = def.indexOf(gene);
  var valueIndex = def.indexOf(value);

  for(let i = 1; i<ar.length-1; i++){
    var res = ar[i].split(split);
    var obj;
    var chrNum = res[chrIndex].replace("chr", "");

    if(localStorage.getItem("genomeType") == "human"){
      if (chrNum == "X"){
        chrNum = "23"
      } else if (chrNum == "Y"){
        chrNum = "24"
      }
    } else if (localStorage.getItem("genomeType") == "mouse"){
      if (chrNum == "X"){
        chrNum = "20"
      } else if (chrNum == "Y"){
        chrNum = "21"
      }
    }
    //entweder pos oder start & end; sample ist optional
    if (posIndex != -1){ //point mutation
      var startPos = res[posIndex];
      if(sampleIndex != -1){ //specified sample
        var samp = res[sampleIndex];
        if (samples.indexOf(samp) == -1){ //sample not in list
          samples.push(samp);
        }
        if (geneIndex != -1) { //specified gene
          var genes = res[geneIndex];

          if(valueIndex != -1){ //specified value
            var values = res[valueIndex];
            obj = {
              sample: samp, chr: chrNum, pos: startPos, gene: genes, value: values
            }
          } else { //no specified value
            obj = {
              sample: samp, chr: chrNum, pos: startPos, gene: genes
            }
          }

        } else { //no specified gene
          if(valueIndex != -1){
            var values = res[valueIndex];
            obj = {
              sample: samp, chr: chrNum, pos: startPos, value: values
            }
          } else {
            obj = {
              sample: samp, chr: chrNum, pos: startPos
            }
          }
        }
      } else { //no specified sample
        if (geneIndex != -1) { //specified gene
          var genes = res[geneIndex];

          if(valueIndex != -1){ //specified value
            var values = res[valueIndex];
            obj = {
              chr: chrNum, pos: startPos, gene: genes, value: values
            }
          } else { //no specified value
            obj = {
              chr: chrNum, pos: startPos, gene: genes
            }
          }
        } else { //no specified gene
          if(valueIndex != -1){ //specified value
            var values = res[valueIndex];
            obj = {
              chr: chrNum, pos: startPos, value: values
            }
          } else { //no specified value
            obj = {
              chr: chrNum, pos: startPos
            }
          }
        }
      }
    } else if (startIndex != -1 && endIndex != -1){ //genomic region
      var startPos = res[startIndex];
      var endPos = res[endIndex];
      if(sampleIndex != -1){ //specified sample
        var samp = res[sampleIndex];
        if (samples.indexOf(samp) == -1){ //sample not in list
          samples.push(samp);
        }
        if (geneIndex != -1){ //specified gene
          var genes = res[geneIndex];

          if(valueIndex != -1){ //specified value
            var values = res[valueIndex];
            obj = {
              sample: samp, chr: chrNum, start: startPos, end: endPos, gene: genes, value: values
            }
          } else { //no specified value
            obj = {
              sample: samp, chr: chrNum, start: startPos, end: endPos, gene: genes
            }
          }
        } else { //no specified gene
          if(valueIndex != -1){ //specified value
            var values = res[valueIndex];
            obj = {
              sample: samp, chr: chrNum, start: startPos, end: endPos, value: values
            }
          } else { //no specified value
            obj = {
              sample: samp, chr: chrNum, start: startPos, end: endPos
            }
          }
        }
      } else { //no specified sample
        if (geneIndex != -1){ //specified gene
          var genes = res[geneIndex];

          if(valueIndex != -1){ //specified value
            var values = res[valueIndex];
            obj = {
              chr: chrNum, start: startPos, end: endPos, gene: genes, value: values
            }
          } else { //no specified value
            obj = {
              chr: chrNum, start: startPos, end: endPos, gene: genes
            }
          }
        } else { //no specified gene
          if(valueIndex != -1){ //specified value
            var values = res[valueIndex];
            obj = {
              chr: chrNum, start: startPos, end: endPos, value: values
            }
          } else { //no specified value
            obj = {
              chr: chrNum, start: startPos, end: endPos
            }
          }
        }
      }
    }
    dataset.push(obj);
  }
  return [dataset, samples];
}

function setConfigurations(data){
  var text = data.split("\n");
  for (let i = 0; i < text.length; i++) {
    const element = text[i];
    var tmp = element.split(" ");
    if (tmp.length > 1){
      localStorage.setItem(tmp[0].replace(":",""), tmp[1]);
    }
  }
}

//--------- CREATE CHECKBOXES OF GIVEN DATASET AND EXPERIMENTS ---------//
function createCheckboxesDataVis(){
  var sample1 = JSON.parse(localStorage.getItem('samples1'));
  var sample2 = JSON.parse(localStorage.getItem('samples2'));
  var data1 = JSON.parse(localStorage.getItem('data1'));
  var data2 = JSON.parse(localStorage.getItem('data2'));

  if (data1.length != 0){
    if (sample1.length == 0){
      createCheckbox("Dataset1", color1[0]);
    } else {
      for(i = 0; i<sample1.length; i++){
        createCheckbox(sample1[i], color1[i]);
      }
    }
  }

  if(data2.length != 0){
    if (sample2.length == 0){
      createCheckbox("Dataset2", color2[0]);
    } else {
      for(i = 0; i<sample2.length; i++){
        createCheckbox(sample2[i], color2[i]);
      }
    }
  }
}


function createCheckbox(id, color){
  var container = document.getElementById('checkboxes')
  var p = document.createElement('p');
  var label = document.createElement('label');
  label.for = id;

  var checkbox = document.createElement('input');
  checkbox.type = "checkbox";
  checkbox.className = "filled-in";
  checkbox.id = id;

  var span = document.createElement('span');
  span.appendChild(document.createTextNode(id));

  var col = document.createElement('span');
  col.style.color = color;
  col.style.fontWeight = "bold";
  col.append(document.createTextNode(' *'));


  label.appendChild(checkbox);
  label.appendChild(span);
  label.appendChild(col);

  p.appendChild(label);

  container.appendChild(p);
}


//--------- CREATE SELECT OF GIVEN DATASET AND EXPERIMENTS ---------//
function createSelectAll(id, selected){
  var sample1 = JSON.parse(localStorage.getItem('samples1'));
  var sample2 = JSON.parse(localStorage.getItem('samples2'));
  var data1 = JSON.parse(localStorage.getItem('data1'));
  var data2 = JSON.parse(localStorage.getItem('data2'));
  if(selected == null)
    document.querySelector("#"+id+" option[value=\"\"]").setAttribute('selected', 'selected');
  if(data1.length != 0)
    createSelect("Dataset1", sample1, id, selected);
  if(data2.length != 0)
    createSelect("Dataset2", sample2, id, selected);
}


function createSelect(dataset, sample, id, selected){
  var container = document.getElementById(id)
  var optgroup = document.createElement('optgroup');
  optgroup.label = dataset;
  if(sample.length == 0){
    var option = document.createElement('option');
    option.value = dataset;
    if(selected != null && selected.includes(dataset)){
      option.selected = true;
    }
    option.appendChild(document.createTextNode(dataset));

    optgroup.appendChild(option);
  } else {
    for (let i=0; i<sample.length; i++){
      var option = document.createElement('option');
      option.value = sample[i];
      if(selected != null && selected.includes(sample[i])){
        option.selected = true;
      }
      option.appendChild(document.createTextNode(sample[i]));

      optgroup.appendChild(option);
    }
  }
  container.appendChild(optgroup);
}


