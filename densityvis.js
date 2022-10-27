
/*
Settings stored in sessionStorage:
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

/*
SLIDER IDs:
- binSize
- bandwith
- amplitude
- numMut
- cluster
- minClusterSize
*/


//fertig
function comparisonDensity(){
  d3.select("#datacomp").selectAll("g").remove();
  d3.select("#datacomp").selectAll("rect").remove();
  d3.select("#datacomp").selectAll("path").remove();
  d3.select("#datacomp").select("text").remove();


  var data1 = JSON.parse(sessionStorage.getItem("data1"));
  var data2 = JSON.parse(sessionStorage.getItem("data2"));

  var genome = sessionStorage.getItem("genomeType"), 
    data_genome,
    centromere;
  
  if(genome == null)
    genome = "human"

  if(genome == "human"){
    data_genome = human_genome;
    centromere = human_centromere;
  }
  else if(genome == "mouse"){
    data_genome = mouse_genome;
    centromere = mouse_centromere;
  }

  //get all selected data --> number of all Mutations
  const dataA = getData(data1, data2,$('#select1').val()), //red
    dataB = getData(data1, data2,$('#select2').val()), //green
    numberOfAllMutationsA = dataA.length,
    numberOfAllMutationsB = dataB.length;
  
  var binSize = sessionStorage.getItem("binSize2"),
    bandwith = sessionStorage.getItem("bandwith2"),
    amplitude = sessionStorage.getItem("amplitude2");

  //get all bins for all chromosomes, get kde-value for each bin 
  var tmpA = getBinDensityAll(dataA, data_genome, numberOfAllMutationsA,$("#datacomp").width(), parseInt(binSize), parseInt(bandwith), parseInt(amplitude)),
    tmpB = getBinDensityAll(dataB,data_genome,numberOfAllMutationsB,$("#datacomp").width(), parseInt(binSize), parseInt(bandwith), parseInt(amplitude));
  var allDensityA = tmpA[1],
    allDensityB = tmpB[1];

  draw_comparison_all(data_genome, centromere, allDensityA, allDensityB)
}


//fertig
function rgbValue(dist1, dist2) {
    var red = 255,
    green = 255;
    if (dist1 == 0 && dist2 == 0){
      return "rgb("+ String(red) + "," + String(green) + ",0";
    } else if (dist1 == 0){
      red = 0;
      green = 255;
    } else if (dist2 == 0) { 
      red = 255;
      green = 0;
    } else {
      if (dist1 > dist2){
        red = 255;
        green = 255/(dist1/dist2);
      }
      else if (dist1 < dist2){
        red = 255/(dist2/dist1);
        green = 255;
      }
    }
    return "rgb("+ String(red) + "," + String(green) + ",0";
}


//fertig
function draw_comparison_distribution(pointsA, pointsB, g, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y){
  /*
    draw the distribution resulting from kde on all chromosomes
  */
  for (let i = 0; i < pointsA.length; ++i) { //points is an array with all bin boundaries stored
    let chrA = pointsA[i] //all bins on chromosome i
    let chrB = pointsB[i] //all bins on chromosome i
    for (let j = 0; j < chrA.length-1   ; ++j) {
      const densityA = (chrA[j][1] + chrA[j+1][1])/2,
        densityB = (chrB[j][1] + chrB[j+1][1])/2,
        opacity = densityA + densityB
      
      /*set polygon, defining the bin*/
      let polygon = getPolygoneOnChromosome(chrA, i, j, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y);
      
      /*draw bin*/
      const bin = g.append("polygon")
        .attr("points", polygon)  
        .style("fill", rgbValue(densityA,densityB))
        .style("opacity", opacity);
    }
  }
}


//fertig
function draw_comparison_all(data_genome, centromere,dataA, dataB)
{
  const width = $("#datacomp").width(),
    height = width*1.2,
    c_width = width/50,
    distance_c = c_width*2,
    c_factor = 320000,
    centromere_x = c_width/5,
    centromere_y = c_width/2;

  const svg = d3
      .select("#datacomp")
      .attr("viewBox", [0, 0,width, height]),
    g = svg.append("g"),
    zoom = d3.zoom()
      .scaleExtent([1, 10])
      .on("zoom", zoomed);

  /*Draw chromosomes depending on genome type*/
  let centromere_d = draw_chromosome(g, data_genome, centromere, distance_c, c_width, c_factor, centromere_x, centromere_y);

  draw_comparison_distribution(dataA, dataB, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y)

  svg.call(zoom);
  
  function zoomed() {
    const {transform} = d3.event;
    g.attr("transform", transform);
    g.attr("stroke-width", 1 / transform.k);
  }
  return svg.node();
}


//fertig
function draw_distribution(dens, numMutations, allGenes, g, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y) {
  /*
    draw the distribution resulting from kde on all chromosomes
  */
  for (let i = 0; i < dens.length; ++i) { //points is an array with all bin boundaries stored
    let chrDens = dens[i] //all bins on chromosome i
    for (let j = 0; j < chrDens.length-1   ; ++j) {
      const opacity = (chrDens[j][1] + chrDens[j+1][1])/2 //set "color" of the bin
      
      /*set polygon, defining the bin*/
      let polygon = getPolygoneOnChromosome(chrDens, i, j, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y);
      
      /*draw grayscaled bin*/
      const bin = g.append("polygon")
        .attr("points", polygon)  
        .on("mouseover", function(d){ //mark bin blue when moving mouse over bin
          d3.select(this)
            .style("fill", "blue")
            .style("opacity", 0.5)
        }) 
        .on("mouseout", function(d) {
          d3.select(this) 
            .transition().duration(100)
            .style("fill", "black")
            .style("opacity", opacity)
        })
        .on("click", function(){ //click on chromosome to change diagram to the chromosome clicked
          sessionStorage.setItem("chrNumber", i)
          all_density()
        })
        .style("fill", "black")
        .style("opacity", opacity);

      bin.append("title") //show detail information of bin, when moving mouse over it
        .text(generate_text(numMutations, allGenes, i, j));
    }
  }
}

//fertig
function getPolygoneOnChromosome(chr, i, j, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y){
  let polygon = "",
    top = (chr[j][0])/(c_factor), //top of the chromosome
    bottom = (chr[j+1][0])/(c_factor), //bottom of the chromosome
    distance_top = 0, //stores difference to normal position, when top of bin in centromere area
    distance_bottom = 0, //stores difference to normal position, when bottom of bin in centromere area
    top_in = false, //is top of bin in centomere area?
    bottom_in = false, //is bottom of bin in centromere area?
    factor = 1; 

  //parts of the bucket are in the area of the centromere
  //top inside 
  if ((top > (centromere[i] - centromere_y)) && (top < (centromere[i] + centromere_y))){ 
    top_in = true;
    factor = -1;
    distance_top = length_mutation_centromere(centromere[i], centromere_x, centromere_y, top);
  }
  //bottom inside
  if ((bottom > (centromere[i] - centromere_y)) && (bottom < (centromere[i] + centromere_y))){
    bottom_in = true;
    distance_bottom = length_mutation_centromere(centromere[i], centromere_x, centromere_y, bottom);
  }
  //(top outside, bottom inside) or (top inside, bottom outside)
  if ((top_in == false && bottom_in == true) || (top_in == true && bottom_in == false)){
    //top outside, bottom inside beneath centromere
    if (top < centromere[i] && bottom > centromere[i] && top_in == false){ //
      polygon = String(i*distance_c + distance_top) + " " + String(top) + ", " + 
        String(i*distance_c) + " " + String(centromere[i] - factor * centromere_y) + ", " +
        String(i*distance_c + centromere_x) + " " + String(centromere[i]) + ", " + 
        String(i*distance_c + distance_bottom) + " " + String(bottom) + ", " + 
        String(i*distance_c + c_width - distance_bottom) + " " + String(bottom) + ", " +
        String(i*distance_c + c_width - centromere_x) + " " + String(centromere[i]) + ", " + 
        String(i*distance_c + c_width) + " " + String(centromere[i] - factor * centromere_y) + ", " +
        String(i*distance_c + c_width - distance_top) + " " + String(top);
      //top inside above centromere, bottom outside
    } else if (top < centromere[i] && bottom > centromere[i] && top_in == true) {
      polygon = String(i*distance_c + distance_top) + " " + String(top) + ", " + 
        String(i*distance_c + centromere_x) + " " + String(centromere[i]) + ", " + 
        String(i*distance_c) + " " + String(centromere[i] - factor * centromere_y) + ", " +
        String(i*distance_c + distance_bottom) + " " + String(bottom) + ", " + 
        String(i*distance_c + c_width - distance_bottom) + " " + String(bottom) + ", " +
        String(i*distance_c + c_width) + " " + String(centromere[i] - factor * centromere_y) + ", " +
        String(i*distance_c + c_width - centromere_x) + " " + String(centromere[i]) + ", " + 
        String(i*distance_c + c_width - distance_top) + " " + String(top);
      //(top ouside and bottom inside above centromere) or (top inside beneath centromere and bottom outside)
    } else {
      polygon = String(i*distance_c + distance_top) + " " + String(top) + ", " + 
        String(i*distance_c) + " " + String(centromere[i] - factor * centromere_y) + ", " +
        String(i*distance_c + distance_bottom) + " " + String(bottom) + ", " + 
        String(i*distance_c + c_width - distance_bottom) + " " + String(bottom) + ", " +
        String(i*distance_c + c_width) + " " + String(centromere[i] - factor * centromere_y) + ", " +
        String(i*distance_c + c_width - distance_top) + " " + String(top);
    }
    // top and bottom inside, centromere between  
  } else if (top_in == true && bottom_in == true && top < centromere[i] && bottom > centromere[i]){
    polygon = String(i*distance_c + distance_top) + " " + String(top) + ", " + 
      String(i*distance_c + centromere_x) + " " + String(centromere[i]) + ", " + 
      String(i*distance_c + distance_bottom) + " " + String(bottom) + ", " + 
      String(i*distance_c + c_width - distance_bottom) + " " + String(bottom) + ", " +
      String(i*distance_c + c_width - centromere_x) + " " + String(centromere[i]) + ", " + 
      String(i*distance_c + c_width - distance_top) + " " + String(top);
    //top and bottom outside, centromere between
  } else if (top_in == false && bottom_in == false && top < centromere[i] && bottom > centromere[i]){
    polygon = String(i*distance_c + distance_top) + " " + String(top) + ", " + 
      String(i*distance_c) + " " + String(centromere[i] - centromere_y) + ", " +
      String(i*distance_c + centromere_x) + " " + String(centromere[i]) + ", " + 
      String(i*distance_c) + " " + String(centromere[i] + centromere_y) + ", " +
      String(i*distance_c + distance_bottom) + " " + String(bottom) + ", " + 
      String(i*distance_c + c_width - distance_bottom) + " " + String(bottom) + ", " +
      String(i*distance_c + c_width) + " " + String(centromere[i] + centromere_y) + ", " +
      String(i*distance_c + c_width - centromere_x) + " " + String(centromere[i]) + ", " + 
      String(i*distance_c + c_width) + " " + String(centromere[i] - centromere_y) + ", " +
      String(i*distance_c + c_width - distance_top) + " " + String(top);
    //top and bottom of bin in the same area (not in centromere area, in centromere area (both above/below centromere))
  } else {
    polygon = String(i*distance_c + distance_top) + " " + String(top) + ", " + 
      String(i*distance_c + distance_bottom) + " " + String(bottom) + ", " + 
      String(i*distance_c + c_width - distance_bottom) + " " + String(bottom) + ", " +
      String(i*distance_c + c_width - distance_top) + " " + String(top);
  }
  return polygon;
}

//points = density --> fertig chrNum startend mit 0
function define_cluster(clusterValue, minClusterSize, points, mutation_bins_gene, mutations_bins_num) {
  /*
    define a region as a cluster, if the kde-value is higher than the threshold value 
    return: array with all clusters, array with all genes effected in one cluster, array with number of mutations per cluster
  */

  var clusterAll = [], //number, chrNum, start, end, mutationNum, gene
    clusterNumber = 1;

  
  for (let chr = 0; chr < points.length; ++chr) { //look at all chromosomes seperately

    let cluster_start = -1,
      cluster_gene = [], //store all effected genes from one cluster
      cluster_mutnumber = 0; //store the number of mutations occured in one cluster
    
    for (let bucket = 0; bucket < points[chr].length-1; ++bucket) { //look at each bucket on chromosome
      if ((points[chr][bucket][1] + points[chr][bucket + 1][1])/2 >= clusterValue){ //kde value higher than threshold 
        cluster_gene.push(mutation_bins_gene[chr][bucket])
        cluster_mutnumber += mutations_bins_num[chr][bucket]
        if (cluster_start == -1){ //firs bucket within a cluster 
          cluster_start = points[chr][bucket][0] //store startpoint
        } else if (bucket == points[chr].length-2){ //last bucket on chromosome --> end of cluster
          let cluster_end = points[chr][bucket][0]; //store endpoint
          if (cluster_end-cluster_start >= minClusterSize){ //everything stored when cluster has minimum size
            let cluster = [clusterNumber, chr, cluster_start, cluster_end, cluster_mutnumber, cluster_gene];
            clusterAll.push(cluster);
          }
          cluster_gene = [];
          cluster_mutnumber = 0;
          cluster_start = -1;
          clusterNumber++;
        }
      } else {
        if (cluster_start != -1) { //there must be a cluster above that bucket (first bucket outside the cluster)
          let cluster_end = points[chr][bucket][0]; 
          if (cluster_end-cluster_start >= minClusterSize){ //everything stored when cluster has minimum size
            let cluster = [clusterNumber, chr, cluster_start, cluster_end, cluster_mutnumber, cluster_gene];
            clusterAll.push(cluster);
          }
          cluster_gene = [];
          cluster_mutnumber = 0;
          cluster_start = -1;
          clusterNumber++;
        }
      }
    }
  }
  return clusterAll
}


//fertig
function draw_cluster(clusterAll, g, distance_c, c_width, c_factor, genomeType) {
  /*
    mark all cluster with an orange box around the bins
    return: cluster number, chromosome number, effected genes, start and end point of each cluster
  */
  var clusterText = "cluster id \tchr num \tstart \tend \tnum of mutations \t genes \tlink to genome browser \n";

  for (let i = 0; i < clusterAll.length; i++) {
    const cluster = clusterAll[i],
      number = cluster[0],
      chrNum = cluster[1],
      start = cluster[2],
      end = cluster[3],
      numMut = cluster[4],
      genes = cluster[5];

    const singleCluster = g.append("polygon") //mark cluster with orange box
      .attr("points", String(chrNum*distance_c) + " " + String(start/c_factor) + ", " +
            String(chrNum*distance_c) + " " + String(end/c_factor) + ", " +
            String(chrNum*distance_c + c_width) + " " + String(end/c_factor) + ", " +
            String(chrNum*distance_c + c_width) + " " + String(start/c_factor))
      .style("stroke", "orange")
      .style("stroke-width", 1.25)
      .style("fill", "none")

    let chrNumber = String(chrNum + 1)
    if(genomeType == "human"){
      if(chrNumber == "23")
        chrNumber = "X"
      else if (chrNumber == "24")
        chrNumber = "Y"
    }
    else if(genomeType == "mouse"){
      if(chrNumber == "20")
        chrNumber = "X"
      else if (chrNumber == "21")
        chrNumber = "Y"
    }

    //lable cluster with an id
    const clusterNumber = g.append("text")
      .attr("x", chrNum*distance_c + c_width*1.5)             
      .attr("y", (start + end)/(2*c_factor) + 5)
      .attr("text-anchor", "middle")  
      .style("font-size", "10px") 
      .text(String(number))

    //set link for the genome browser
    var link = linkGenomeBrowser(chrNumber, start, end);

    //link cluster to the UCSC genome browser (click on cluster id)
    clusterNumber.html('<a href= "'+link+'" target="_blank">' + String(number) + '</a>')


    //get all effected gene in that cluster and sort them
    var cluster_gene = []
    for (let k = 0; k < genes.length; ++k) {
      for (let l = 0; l < genes[k].length; ++l) {
        cluster_gene.push(genes[k][l])
      }
    }
    cluster_gene = cluster_gene.filter(function(ele , pos){
      return cluster_gene.indexOf(ele) == pos;
    }) 
    cluster_gene.sort()

    //show detail information of cluster, when moving mouse over id
    if (cluster_gene.length > 0){ //cluster has effected genes
      //chromosome number, position on chromosome, number of mutations, effected genes
      clusterNumber.append("title") 
        .text("chromosome: " + chrNumber + "\nposition: " + String(start) + "-" + String(end) + 
              "\nnumber of mutations: " + String(numMut) + "\ngene: " + String(cluster_gene));
      clusterText += String(number) + "\t" + chrNumber + "\t" + String(start) + "\t" + String(end) + "\t" + String(numMut) + "\t" + String(cluster_gene) + "\t" + String(link) + "\n"; 
    } else { //cluster has no effected genes
      clusterNumber.append("title") 
        .text("chromosome: " + chrNumber + "\nposition: " + String(start) + "-" + String(end) + 
              "\nnumber of mutations: " + String(numMut) + "\nno specific gene");
      clusterText += String(number) + "\t" + chrNumber + "\t" + String(start) + "\t" + String(end) + "\t" + String(numMut) + "\t" + "no sepecific" + "\t" + String(link) + "\n";
    }

    
  }
  sessionStorage.setItem("clusterText", clusterText);
}


//fertig
function mark_bins(bins, mutationNumbers, numMutationsThreshold, g, distance_c, c_width, c_factor) {
  //mutationNumbers is an array with the number of mutations in each bin
  //numMutations is the number of mutations a bin must have to be marked

  /*mark bins with red line on right side of the bin, when bin has more than numMutations mutations*/
  for (let i = 0; i < bins.length; ++i) {
    for (let j = 0; j < bins[i].length; ++j) {
      if (mutationNumbers[i][j] >= numMutationsThreshold){ //mark the bin
        g.append("polygon")
          .attr("points", String(i*distance_c + c_width + 2) + " " + String(bins[i][j][0]/c_factor) + ", " +
                String(i*distance_c + c_width + 2) + " " + String(bins[i][j+1][0]/c_factor))
          .style("stroke", "red")
          .style("stroke-width", 2)
      }
    }
  }
}

//fertig
function linkGenomeBrowser(chrNumber, start, end){ 
  /*
    generates the link zu the ucsc genome browser
    returns the url as a string
  */

  return "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=" + sessionStorage.getItem("genomeRef") + "&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr" + String(chrNumber) + "%3A" + String(start) + "%2D" + String(end) + "&hgsid=287468263_f4ah5As0fHokqKgjwt2PaNwTeBNM"  
}


//fertig
function generate_text(numMutations, allGenes, i, j) {
  /*
    generates the text for each bin, shown when moving mouse over bin
  */
  if (allGenes[i][j].length > 0){
    return "mutations: "  + String(numMutations[i][j]) + "\n"  + "gene: " + allGenes[i][j]
  } else {
    return "mutations: "  + String(numMutations[i][j]) 
  }
}



//fertig
function draw_density_all(data_genome, centromere, data, density, numMutations, allGenes, cluster, genomeType, minNumMut) {
  const width = $("#densvisAll").width(),
    height = width*1.2,
    c_width = width/50,
    distance_c = c_width*2,
    c_factor = 320000,
    centromere_x = c_width/5,
    centromere_y = c_width/2;

  const svg = d3
      .select("#densvisAll")
      .attr("viewBox", [0, 0,width, height]),
    g = svg.append("g"),
    zoom = d3.zoom()
      .scaleExtent([1, 10])
      .on("zoom", zoomed);
  

  /*Draw chromosomes depending on genome type*/
  let centromere_d = draw_chromosome(g, data_genome, centromere, distance_c, c_width, c_factor, centromere_x, centromere_y);
  
  /*draw distribution*/
  if (document.getElementById('distributionCbx').checked){
    draw_distribution(density, numMutations, allGenes, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y)
  }

  /*draw mutation*/
  if (document.getElementById('mutationCbx').checked){ 
    draw_all_mutation_points(data, g, distance_c, c_width, c_factor, data_genome)
  }
  
  /*mark bins*/
  if (document.getElementById('binMutCbx').checked){
    mark_bins(density, numMutations, minNumMut, g, distance_c, c_width, c_factor);
  }

  /*mark cluster*/
  if (document.getElementById('clusterCbx').checked){
    draw_cluster(cluster, g, distance_c, c_width, c_factor, genomeType)
  }
  
  svg.call(zoom);
  
  function zoomed() {
    const {transform} = d3.event;
    g.attr("transform", transform);
    g.attr("stroke-width", 1 / transform.k);
  }

  return svg.node();
}


//fertig
function draw_diagram_chr(chromosome, data_genome, bins, density) {
  const width = $("#densvisChr").width(), //width of the diagram
    height = width/5; //height of the diagram

  /*data/mutations*/
  const length = data_genome[chromosome-1].size;  //length of chromosome from given number

  const svg = d3
    .select("#densvisChr")
    .attr("width", width)
    .attr("height", height),
    margin = {top: 20, right: 40, bottom: 30, left: 25};

  const x = d3.scaleLinear()
    .domain([0,length])
    .range([margin.left, width - margin.right]);

  var tmp = [];
  for(let i=0; i<density.length; i++){
    tmp.push(density[i][0]);
  }
  
  function color(d){
    if (document.getElementById('clusterCbx').checked){
      var index = tmp.indexOf(d.x0);
      if (index < density.length-2){
        if ((density[index][1] + density[index+1][1])/2 >= $('#cluster').val())
          return "orange"
      }
    } 
    return "#bbb"
  }
  
  const y = d3.scaleLinear()
    .domain([0, d3.max(bins, d => d.length)])
    .range([height - margin.bottom, margin.top]);

  /*set title (chromosome id)*/        
  svg.append("text")
    .attr("x", (width / 2))             
    .attr("y", 5 + margin.top / 2)
    .attr("text-anchor", "middle")  
    .style("font-size", "9x") 
    .text("Chromosome " + String(chromosome));

  
  /*draw coordinatesystem*/
  svg.append("g")
    .attr("class", "axis axis--x")
    .attr("transform", "translate(0," + (height - margin.bottom) + ")")
    .call(d3.axisBottom(x).tickFormat(d3.format(".1e")));
  svg.append("g")
    .attr("class", "axis axis--y")
    .attr("transform", "translate(" + margin.left + ",0)")
    .call(d3.axisLeft(y).ticks(5));


  /*draw histogramm of mutations*/
  svg.selectAll("rect")
    .data(bins)
    .enter()
    .append("rect")
      .attr("x", function(d) { return x(d.x0) + 1; })
      .attr("y", function(d) { return y(d.length); })
      .attr("width", function(d) { return x(d.x1) - x(d.x0); })
      .attr("height", function(d) { return y(0) - y(d.length); })
      .style("fill", function(d) { return color(d); });


  /*draw kde graph*/
  var path = svg.append("path")
      .datum(density)
      .attr("fill", "none")
      .attr("stroke", "#000")
      .attr("stroke-width", 1.5)
      .attr("stroke-linejoin", "round")
      .attr("d",  d3.line() 
          .curve(d3.curveBasis)
          .x(function(d) { return x(d[0])})
          .y(function(d) { return y(d[1])}));

  /*add threshold graph for cluster when checkbox is set*/
  if (document.getElementById('clusterCbx').checked){
    var cluster = parseFloat(($('#cluster').val()))
    const clusterThreshold = [[0,cluster],[length,cluster]]
    var clusterLine = svg.append("path")
      .datum(clusterThreshold)
      .attr("stroke", "orange")
      .attr("stroke-width", 1)
      .attr("stroke-linejoin", "bevel")
      .attr("d",  d3.line()
          .curve(d3.curveBasis)
          .x(d => x(d[0]))
          .y(d => y(d[1])));
  }

  return svg.node();
}


//fertig
function kde(kernel, thresholds, data, numberOfAllMutations) {
  if (data.length == 0)
    return thresholds.map(t => [t,0]);
  return thresholds.map(t => [t, (d3.sum(data, d => kernel(t - d))/numberOfAllMutations)]); //get kde value at the bin boundaries
}

//fertig
function epanechnikov(bw, amplitude) {
  return function(v) {
    return Math.abs(v /= bw) <= 1 ? amplitude * (1 - v * v) / bw : 0;
  };
}

//fertig
function filter_mutations(chr_number, data){ //filter dataset by mutations
  const data_chr = data.filter(x => x.chr === chr_number);
  let res = [];
  for (let i = 0; i < data_chr.length; ++i) {
    var mut = parseInt(data_chr[i].pos);
    if (isNaN(mut))
      mut = parseInt(data_chr[i].start)
    res.push(mut);
  }
  return res
}


//fertig
function all_density(){
  d3.select("#densvisChr").selectAll("g").remove();
  d3.select("#densvisChr").selectAll("rect").remove();
  d3.select("#densvisChr").selectAll("path").remove();
  d3.select("#densvisChr").select("text").remove();

  d3.select("#densvisAll").selectAll("g").remove();
  d3.select("#densvisAll").selectAll("rect").remove();
  d3.select("#densvisAll").selectAll("path").remove();
  d3.select("#densvisAll").select("text").remove();

  var data1 = JSON.parse(sessionStorage.getItem("data1"));
  var data2 = JSON.parse(sessionStorage.getItem("data2"));

  var genome = sessionStorage.getItem("genomeType"), 
    data_genome,
    centromere;

  if (genome == null) //genome wurde noch nicht gesetzt
    genome = "human";

  if(genome == "human"){
    data_genome = human_genome;
    centromere = human_centromere;
  }
  else if(genome == "mouse"){
    data_genome = mouse_genome;
    centromere = mouse_centromere;
  }

  //get all selected data --> number of all Mutations
  const data = getData(data1, data2, $('#select').val()),
    numberOfAllMutations = data.length;

  //get all bins for all chromosomes, get kde-value for each bin 
  //Parameter aus Speicher lesen
  var binSize = sessionStorage.getItem("binSize1"),
    bandwith = sessionStorage.getItem("bandwith1"),
    amplitude = sessionStorage.getItem("amplitude1"),
    clusterThreshold = sessionStorage.getItem("cluster"),
    minClusterSize = sessionStorage.getItem("minClusterSize"),
    minNumMut = sessionStorage.getItem("numMut")
  var tmp = getBinDensityAll(data, data_genome, numberOfAllMutations, $("#densvisChr").width(), parseInt(binSize), parseInt(bandwith), parseInt(amplitude));
  var allBins = tmp[0],
    allDensity = tmp[1],
    allMutNum = tmp[2],
    allGene = tmp[3];

  var numberChr = parseInt(sessionStorage.getItem("chrNumber"))
  if(isNaN(numberChr))
    numberChr = 0

  var cluster = define_cluster(parseFloat(clusterThreshold), parseInt(minClusterSize), allDensity, allGene, allMutNum);

  //draw chromosome diagram of clicked chromosome
  draw_diagram_chr(numberChr+1, data_genome, allBins[numberChr], allDensity[numberChr])

  //var all_chr_diagrams = all_diagramms(data_genome, $('#select').val() ,data1, data2);
  draw_density_all(data_genome, centromere, data, allDensity, allMutNum, allGene, cluster, genome, parseInt(minNumMut))
}


//fertig
function getData(data1, data2, select){ //get data from selected checkboxes
  var sample1 = JSON.parse(sessionStorage.getItem('samples1'));
  var sample2 = JSON.parse(sessionStorage.getItem('samples2'));
   
  var data = [];
  for (i = 0; i<select.length; i++){
    if(select[i] == "Dataset1")
      data = data.concat(data1);
    else if (select[i] == "Dataset2")
      data = data.concat(data2);
    else{
      let ind1 = sample1.indexOf(select[i]);
      let ind2 = sample2.indexOf(select[i])
      if (ind1 != -1){
        data = data.concat(data1.filter(x => x.sample === select[i])); 
      } else if (ind2 != -1){
        data = data.concat(data2.filter(x => x.sample === select[i])); 
      }
    }
  }
  return data;
}

//fertig
function getBinDensity(data, chromosome, data_genome, numberOfAllMutations, width, binSize, bandwidth, amplitude){
  //return all bins mit included elements and retruns density at bin boundaries
  const length = data_genome[chromosome-1].size;  //length of chromosome from given number starting with 0
  let chr = filter_mutations(String(chromosome), data);   //all mutations of the given chromosome

  var margin = {top: 20, right: 85, bottom: 30, left: 25}

  const x = d3.scaleLinear()
  .domain([0,length])
  .range([margin.left, width - margin.right]);

    /*get bin boundaries that have the same distance (exept the last bin), with the same bin size for all chromosomes*/
  function equiSizedBins() {
    var binBounds = [];
    for (let j = 0; j < length; j += binSize) { 
      binBounds.push(j);
    }
    if(binBounds[binBounds.length - 1] != length){ //add last remaining bin (smaller size)
      binBounds.push(parseInt(length))
    }
    return binBounds;
  }

  var thresholds = equiSizedBins(), //array with boundaries of the bins
    bins = d3.histogram().domain(x.domain()).thresholds(thresholds)(chr),
    density = kde(epanechnikov(bandwidth, amplitude), thresholds, chr, numberOfAllMutations); //get density of all mutations using kde
  bins.pop()

  /* save mutation number and gene for each bin */
  var mutationNum = [],
    gene = [];
  for (let i = 0; i < bins.length; i++) {
    mutationNum.push(bins[i].length);

    let geneBin = [];
    for (let j = 0; j < bins[i].length; j++) {
      const mutation = bins[i][j];
      const data_chr = data.filter(x => x.chr === String(chromosome) && x.pos === String(mutation)),
        data_gene = data_chr.map(data_chr => data_chr.gene);
      geneBin.push(data_gene[0])
    }
    //remove duplicates
    const gene_wo_d = geneBin.filter(function(ele , pos){
      return geneBin.indexOf(ele) == pos;
    }) 
    //sort
    gene_wo_d.sort()
    //remove N/A, n/a, ""
    if (gene_wo_d.includes("N/A")){
      gene_wo_d.splice(gene_wo_d.indexOf("N/A"),1)  //remove
    } 
    if (gene_wo_d.includes("n/a")){
      gene_wo_d.splice(gene_wo_d.indexOf("n/a"),1)  //remove
    } 
    if (gene_wo_d.includes("")){
      gene_wo_d.splice(gene_wo_d.indexOf(""),1)  //remove
    }
    if (gene_wo_d.includes(undefined)){
      gene_wo_d.splice(gene_wo_d.indexOf(undefined),1)  //remove
    }
    //add gene from bucket to chromosome
    gene.push(gene_wo_d)
  }

  return [bins, density, mutationNum, gene]; //mutations number = length of each bin
}

//fertig --> Bins, Density, NumMutations, Genes
function getBinDensityAll(data, data_genome, numberOfAllMutations, width, binSize, bandwidth, amplitude){
  var binAll =[],
    densityAll =[],
    mutationNumAll = [],
    geneAll = [];

  for (chrNum = 1; chrNum < data_genome.length+1; chrNum++){
    var tmp = getBinDensity(data, chrNum, data_genome, numberOfAllMutations, width, binSize, bandwidth, amplitude);
    binAll.push(tmp[0]);
    densityAll.push(tmp[1]);
    mutationNumAll.push(tmp[2]);
    geneAll.push(tmp[3]);
  }
  return [binAll, densityAll, mutationNumAll, geneAll]
}
