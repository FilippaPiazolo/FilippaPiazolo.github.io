
/*
Settings stored in localStorage:

Genom Type: #genomeType
Genome Reference: #genomeRef
Dataset 1: #data1 
Dataset 2: #data2
Data Type 1: #datatype1
Data Type 2: #datatype2
Number of mutations dataset 1: #numdata1
Number of mutations dataset 2: #numdata2
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


function _densData(select,getDensityOptions,samples){return(
select({
  title: "Data for distribution",
  description: "Please select the data used for the distribution and clustering.",
  options: getDensityOptions(samples),
  value: "all data"
})
)}

function _visions(checkbox,numMut){return(
checkbox({
  title: "",
  description: "Please select distribution and click on a chromosome to see its diagram and click on checkbox for showing elements on all chromosomes",
  options: [
    { value: "distribution", label: "distribution", color: "black"},
    { value: "mutation", label: "mutation (one color)", color: "red"},
    { value: "bins", label: "show bins with " + String(numMut) + " or more mutations", color: "green"},
    { value: "cluster", label: "show cluster", color: "orange"}
  ],
  value: ["distribution"]
})
)}

function _19(all_diagramms,chrNumber){return(
all_diagramms[chrNumber][1]
)}

function _20(resultsCluster){return(
resultsCluster[0]
)}

function _21(generate_text_cluster,resultsCluster,md)
{
  let text = generate_text_cluster(resultsCluster[1])
  if ( text != ""){
    return md`
### Cluster and Gene
| Cluster (Browser)| Chromosome | Gene  | Position  |
| ------- | ---------- | ----- | -------------- |
${text}
    `
  } else {
    return md`select checkbox *show cluster* to see a list of all cluster`
  }
}


function _densData1(select,getDensityOptions,samples){return(
select({
  title: "Data1",
  description: "marked red",
  options: getDensityOptions(samples),
  value: "all data"
})
)}

function _densData2(select,getDensityOptions,samples){return(
select({
  title: "Data2",
  description: "marked green",
  options: getDensityOptions(samples),
  value: "all data"
})
)}

function _comparison(d3,genome,human_centromere,human_genome,mouse_centromere,mouse_genome,data,draw_diagram,getBinBoarders,densData1,densData2,draw_chromosome,draw_comparison_distribution)
{
  var width = 1000,
    height = 1000,
    c_width = 15,
    distance_c = c_width + 20,
    c_factor = 300000,
    centromere_x = 3,
    centromere_y = 7,
    all_cluster = [];

  const svg = d3
      .create('svg')
      .attr("viewBox", [0, 0,width, height]),
    g = svg.append("g"),
    zoom = d3.zoom()
      .scaleExtent([1, 10])
      .on("zoom", zoomed);

  let centromere, data_genome;
  if (genome == "human"){
    centromere = human_centromere;
    data_genome = human_genome;
  } else if (genome == "mouse"){
    centromere = mouse_centromere;
    data_genome = mouse_genome;
  }

  function getData(densData){
    let tmp = [];
    for (let i = 1; i <= data_genome.length; ++i) {
      let experiment = [];
      if (densData == "all data"){
        for (let j=0; j<data.length; j++){
          experiment = experiment.concat(data[j]);
        }
      } 
      else if (densData == "dataset1")
        experiment = experiment.concat(data[0]);
      else if (densData == "dataset2")
        experiment = experiment.concat(data[1]);
      else{
        for (let j=0; j<data.length; j++){
          experiment = experiment.concat(data[j].filter(x => x.Sample === densData));
        }
      }
      tmp.push(draw_diagram(experiment, i, data_genome));
    }
    //array with svg-elements for diagrams of all chromosomes
    return tmp
  }

  let dataA = getBinBoarders(getData(densData1)),
    dataB = getBinBoarders(getData(densData2));
  

  /*Draw chromosomes depending on genome type*/
  let centromere_d = draw_chromosome(g, data_genome, centromere, distance_c, c_width, c_factor, centromere_x, centromere_y);

  draw_comparison_distribution(dataA, dataB, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y)

  svg.call(zoom);
  
  function zoomed(event) {
    const {transform} = event;
    g.attr("transform", transform);
    g.attr("stroke-width", 1 / transform.k);
  }
  return svg.node();
}







function _chrNumber(){return(
0
)}

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
      
      /*draw grayscaled bin*/
      const bin = g.append("polygon")
        .attr("points", polygon)  
        .style("fill", rgbValue(densityA,densityB))
        .style("opacity", opacity);
    }
  }
}


function draw_distribution(points, mut_bins, g, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y) {
  /*
    draw the distribution resulting from kde on all chromosomes
  */
  for (let i = 0; i < points.length; ++i) { //points is an array with all bin boundaries stored
    let chr = points[i] //all bins on chromosome i
    for (let j = 0; j < chr.length-1   ; ++j) {
      const opacity = (chr[j][1] + chr[j+1][1])/2 //set "color" of the bin
      
      /*set polygon, defining the bin*/
      let polygon = getPolygoneOnChromosome(chr, i, j, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y);
      
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
          $0.value = i
        })
        .style("fill", "black")
        .style("opacity", opacity);

      bin.append("title") //show detail information of bin, when moving mouse over it
        .text(generate_text(mut_bins, i, j));
    }
  }
}


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


function define_cluster(clusterValue, points, mutation_bins_gene, mutations_bins_num) {
  /*
    define a region as a cluster, if the kde-value is higher than the threshold value 
    return: array with all clusters, array with all genes effected in one cluster, array with number of mutations per cluster
  */
  let all_cluster = [], //store all clusters (start and end point)
    all_cluster_gene = [], //store all genes effected in cluster
    all_cluster_mutnumber = [] //store number of mutations per cluster
  
  for (let chr = 0; chr < points.length; ++chr) { //look at all chromosomes seperately
    let chr_cluster = [], //store all cluster per chromosome
      chr_cluster_gene = [], //store all effected genes per cluster per chromosome
      chr_cluster_mutnumber = [], //store number of mutations per cluster per chromosome
      cluster_start = -1,
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
          if (cluster_end-cluster_start > minClusterSize){ //everything stored when cluster has minimum size
            let cluster = [cluster_start, cluster_end] 
            chr_cluster.push(cluster) //store cluster (beginning, end)
            chr_cluster_gene.push(cluster_gene) //store all effected genes from that cluster
            chr_cluster_mutnumber.push(cluster_mutnumber) //store all mutations from that cluster
          }
          cluster_gene = []
          cluster_mutnumber = 0
          cluster_start = -1
        }
      } else {
        if (cluster_start != -1) { //there must be a cluster above that bucket (first bucket outside the cluster)
          let cluster_end = points[chr][bucket][0]; 
          if (cluster_end-cluster_start > minClusterSize){ //everything stored when cluster has minimum size
            let cluster = [cluster_start, cluster_end]
            chr_cluster.push(cluster)
            chr_cluster_gene.push(cluster_gene)
            chr_cluster_mutnumber.push(cluster_mutnumber)
          }
          cluster_gene = []
          cluster_mutnumber = 0
          cluster_start = -1
        }
      }
    }
    all_cluster.push(chr_cluster)
    all_cluster_gene.push(chr_cluster_gene)
    all_cluster_mutnumber.push(chr_cluster_mutnumber)
  }
  return [all_cluster,all_cluster_gene, all_cluster_mutnumber] 
}


function draw_cluster(points, mut_bins, g, distance_c, c_width, c_factor, clusterValue) {
  /*
    mark all cluster with an orange box around the bins
    return: cluster number, chromosome number, effected genes, start and end point of each cluster
  */
  const all = define_cluster(clusterValue, points, mut_bins[1], mut_bins[0]), 
    all_cluster = all[0],
    all_cluster_gene = all[1],
    all_cluster_mutNum = all[2],
    res = []
  var number = 1

  /*draw cluster*/
  for (let i = 0; i < all_cluster.length; ++i) {
    for (let j = 0; j < all_cluster[i].length; ++j) {
      const cluster = g.append("polygon") //mark cluster with orange box
        .attr("points", String(i*distance_c) + " " + String(all_cluster[i][j][0]/c_factor) + ", " +
              String(i*distance_c) + " " + String(all_cluster[i][j][1]/c_factor) + ", " +
              String(i*distance_c + c_width) + " " + String(all_cluster[i][j][1]/c_factor) + ", " +
              String(i*distance_c + c_width) + " " + String(all_cluster[i][j][0]/c_factor))
        .style("stroke", "orange")
        .style("stroke-width", 1)
        .style("fill", "none")
        
      let chrNumber = String(i + 1)
      if(genome == "human"){
        if(chrNumber == "23")
          chrNumber = "X"
        else if (chrNumber == "24")
          chrNumber = "Y"
      }
      else if(genome == "mouse"){
        if(chrNumber == "20")
          chrNumber = "X"
        else if (chrNumber == "21")
          chrNumber = "Y"
      }

      //lable cluster with an id
      const clusterNumber = g.append("text")
        .attr("x", i*distance_c + c_width*1.5)             
        .attr("y", (all_cluster[i][j][0] + all_cluster[i][j][1])/(2*c_factor) + 5)
        .attr("text-anchor", "middle")  
        .style("font-size", "10px") 
        .text(String(number))

      //set link for the genome browser
      var link = linkGenomeBrowser(chrNumber, all_cluster[i][j][0], all_cluster[i][j][1]);

      //link cluster to the UCSC genome browser (click on cluster id)
      clusterNumber.html('<a href= "'+link+'" target="_blank">' + String(number) + '</a>')


      //get all effected gene in that cluster and sort them
      var cluster_gene = []
      for (let k = 0; k < all_cluster_gene[i][j].length; ++k) {
        for (let l = 0; l < all_cluster_gene[i][j][k].length; ++l) {
          cluster_gene.push(all_cluster_gene[i][j][k][l])
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
          .text("chr" + chrNumber + ";" + String(all_cluster[i][j][0]) + "-" + String(all_cluster[i][j][1]) + 
                "\n number of mutations: " + String(all_cluster_mutNum[i][j]) + "\n gene: " + String(cluster_gene)) 
      } else { //cluster has no effected genes
        clusterNumber.append("title") 
          .text("chr" + chrNumber + ";" + String(all_cluster[i][j][0]) + "-" + String(all_cluster[i][j][1]) + 
                "\n number of mutations: " + String(all_cluster_mutNum[i][j]) + "\n no specific gene")
      }
      
      res.push([number, chrNumber, cluster_gene, [all_cluster[i][j][0], all_cluster[i][j][1]]]) 
      number++
    }
  }
  return res
}


function mark_bins(points, mutationNumbers, numMutations, g, distance_c, c_width, c_factor) {
  //mutationNumbers is an array with the number of mutations in each bin
  //numMutations is the number of mutations a bin must have to be marked

  /*mark bins with red line on right side of the bin, when bin has more than numMutations mutations*/
  for (let i = 0; i < points.length; ++i) {
    for (let j = 0; j < points[i].length; ++j) {
      if (mutationNumbers[i][j] >= numMutations){ //mark the bin
        g.append("polygon")
          .attr("points", String(i*distance_c + c_width + 2) + " " + String(points[i][j][0]/c_factor) + ", " +
                String(i*distance_c + c_width + 2) + " " + String(points[i][j+1][0]/c_factor))
          .style("stroke", "green")
          .style("stroke-width", 2)
      }
    }
  }
}


function linkGenomeBrowser(chrNumber, start, end){ 
  /*
    generates the link zu the ucsc genome browser
    returns the url as a string
  */
  var setRef = "";
  if(refGenome == "GRCh37/hg19")
    setRef = "hg19";
  else if(refGenome == "GRCh38/hg38")
    setRef = "hg38";
  else if(refGenome == "GRCm38/mm10")
    setRef = "mm10";
  else if(refGenome == "GRCm39/mm39")
    setRef = "mm39";

  return "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=" + setRef + "&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr" + String(chrNumber) + "%3A" + String(start) + "%2D" + String(end) + "&hgsid=287468263_f4ah5As0fHokqKgjwt2PaNwTeBNM"  
}


function generate_text(mut_bins, i, j) {
  /*
    generates the text for each bin, shown when moving mouse over bin
  */
  if (mut_bins[1][i][j].length > 0){
    return "mutations: "  + String(mut_bins[0][i][j]) + "\n"  + "gene: " + mut_bins[1][i][j]
  } else {
    return "mutations: "  + String(mut_bins[0][i][j]) 
  }
}


function generate_text_cluster(allCluster) { 
  /*
    generates text for table with overview of resulting clusters
  */
  var res = ""
  for (let cluster = 0; cluster < allCluster.length; ++cluster) {
    var gene = ""
    if (allCluster[cluster][2].length > 0)
      gene = String(allCluster[cluster][2])
    else
      gene = "no specific gene"
    res += "|" + "[" + String(allCluster[cluster][0])+"]("+ linkGenomeBrowser(allCluster[cluster][1], allCluster[cluster][3][0], allCluster[cluster][3][1]) + ")" + "|" + String(allCluster[cluster][1]) + 
      "|" + gene + "|" + String(allCluster[cluster][3][0]) + "-" + String(allCluster[cluster][3][1]) + "|\n"
  }
  return res
}


function getBinBoarders(diagramms){
  let tmp = []
  for (let i = 0; i < diagramms.length; ++i) {
    tmp.push(diagramms[i][0])
  }
  return tmp
}


function getMutationsGeneNumber(data, diagramms) {
  let all = [],
    gene_all = []
  for (let chrNum = 0; chrNum < diagramms.length; ++chrNum) { //all chromosomes
    let chr = [],
      gene_chr = []
    for (let j = 0; j < diagramms[chrNum][2].length; ++j) { //all buckets on each chromosome
      chr.push(diagramms[chrNum][2][j].length) //number of mutations in each bucket on each chromosome
      //specific genes in each bucket
      let gene = []
      for (let k = 0; k < diagramms[chrNum][2][j].length; ++k) {
        let mutation = diagramms[chrNum][2][j][k]
        const data_chr = data.filter(x => x.chr === String(chrNum+1) && x.pos === String(mutation)),
          data_gene = data_chr.map(data_chr => data_chr.GENE)
        gene.push(data_gene[0])
      }
      //remove duplicates
      const gene_wo_d = gene.filter(function(ele , pos){
        return gene.indexOf(ele) == pos;
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
      gene_chr.push(gene_wo_d)
    }
    all.push(chr)
    gene_all.push(gene_chr)
  }
  return [all, gene_all]
}


function all_together(data, diagramms) {
  var width = 1000,
    height = 1000,
    c_width = 15,
    distance_c = c_width + 20,
    c_factor = 300000,
    centromere_x = 3,
    centromere_y = 7,
    all_cluster = [];

  const svg = d3
      .create('svg')
      .attr("viewBox", [0, 0,width, height]),
    g = svg.append("g"),
    zoom = d3.zoom()
      .scaleExtent([1, 10])
      .on("zoom", zoomed);

  let centromere, data_genome;
  if (genome == "human"){
    centromere = human_centromere;
    data_genome = human_genome;
  } else if (genome == "mouse"){
    centromere = mouse_centromere;
    data_genome = mouse_genome;
  }

  let allBoarders = getBinBoarders(all_diagramms),
    mutation_gene_num = getMutationsGeneNumber(data, diagramms);
  

  /*Draw chromosomes depending on genome type*/
  let centromere_d = draw_chromosome(g, data_genome, centromere, distance_c, c_width, c_factor, centromere_x, centromere_y);
  
  /*draw distribution*/
  if (visions.includes("distribution")){
    draw_distribution(allBoarders, mutation_gene_num, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y)
  }
  /*draw mutation*/
  if (visions.includes("mutation")){ 
    draw_all_mutation_points(data, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y, data_genome)
  }
  /*draw colored mutations*/
 /* if (visions.includes("mutation2")){
    draw_mutations(data, color, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y, samples, data_genome);
  }*/
  /*mark bins*/
  if (visions.includes("bins")){
    mark_bins(allBoarders, mutation_gene_num[0], numMut, g, distance_c, c_width, c_factor);
  }
  /*mark cluster*/
  if (visions.includes("cluster")){
    all_cluster = draw_cluster(allBoarders, mutation_gene_num, g, distance_c, c_width, c_factor, cluster)
  }
  
  svg.call(zoom);
  
  function zoomed(event) {
    const {transform} = event;
    g.attr("transform", transform);
    g.attr("stroke-width", 1 / transform.k);
  }
  return [svg.node(), all_cluster];
}


function _resultsCluster(densData,data,all_together,all_diagramms)
{
  if (densData == "all data"){
    let experiment = [];
    for (let i=0; i<data.length;i++){
      experiment = experiment.concat(data[i])
    }
    return all_together(experiment, all_diagramms)
  }
  else if (densData == "dataset1")
    return all_together(data[0], all_diagramms)
  else if (densData == "dataset2")
    return all_together(data[1], all_diagramms)
  else{
    const experiment = data.filter(x => x.Sample === densData); 
    return all_together(experiment, all_diagramms)
  }
}



function _all_diagramms(genome,human_genome,mouse_genome,densData,data,draw_diagram)
{
  let tmp = [],
    data_genome;
  if(genome == "human")
    data_genome = human_genome;
  if(genome == "mouse")
    data_genome = mouse_genome;
  for (let i = 1; i <= data_genome.length; ++i) {
    let experiment = [];
    if (densData == "all data"){
      for (let j=0; j<data.length; j++){
        experiment = experiment.concat(data[j]);
      }
    } 
    else if (densData == "dataset1")
      experiment = experiment.concat(data[0]);
    else if (densData == "dataset2")
      experiment = experiment.concat(data[1]);
    else{
      for (let j=0; j<data.length; j++){
        experiment = experiment.concat(data[j].filter(x => x.Sample === densData));
      }
    }
    tmp.push(draw_diagram(experiment, i, data_genome));
  }
  //array with svg-elements for diagrams of all chromosomes
  return tmp
}


function draw_diagram(data, chromosome, data_genome) {
  const width = 900, //width of the diagram
    height = 170; //height of the diagram

  /*data/mutations*/
  const length = data_genome[chromosome-1].size;  //length of chromosome from given number
  let chr = filter_mutations(String(chromosome), data),   //all mutations of the given chromosome
    n = chr.length;   //number of mutations

  const svg = d3
    .create("svg")
    .attr("width", width)
    .attr("height", height),
    margin = {top: 20, right: 85, bottom: 30, left: 40};

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

  const thresholds = equiSizedBins(), //array with boundaries of the bins
    bins = d3.histogram().domain(x.domain()).thresholds(thresholds)(chr),
    density = kde(epanechnikov(bandwidth), thresholds, chr); //get density of all mutations using kde
  //console.log(density);
  var tmp = [];
  for(let i=0; i<density.length; i++){
    tmp.push(density[i][0]);
  }
  
  function color(d){
    if (visions.includes("cluster")){
      var index = tmp.indexOf(d.x0);
      if (index < density.length-2){
        if ((density[index][1] + density[index+1][1])/2 >= cluster)
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
    .call(d3.axisBottom(x));
  svg.append("g")
    .attr("class", "axis axis--y")
    .attr("transform", "translate(" + margin.left + ",0)")
    .call(d3.axisLeft(y)
    .ticks(5));


  /*draw histogramm of mutations*/
  svg.insert("g", "*")
    .selectAll("rect")
    .data(bins)
    .join(enter => enter.append("rect")
      .attr("x", function(d) {
        return x(d.x0) + 1; })
      .attr("y", function(d) { return y(d.length); })
      .attr("width", function(d) { return x(d.x1) - x(d.x0) - 1; })
      .attr("height", function(d) { return y(0) - y(d.length); })
      .style("fill", function(d) { return color(d); }));


  /*draw kde graph*/
  var path = svg.append("path")
      .datum(density)
      .attr("fill", "none")
      .attr("stroke", "#000")
      .attr("stroke-width", 1.5)
      .attr("stroke-linejoin", "round")
      .attr("d",  d3.line()
          .curve(d3.curveBasis)
          .x(d => x(d[0]))
          .y(d => y(d[1])));

  /*add threshold graph for cluster when checkbox is set*/
  if (visions.includes("cluster")){
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

  //return kde-value at the boarder of each bin, all diagramms, all bins
  return [path._groups[0][0].__data__, svg.node(), bins];
}


function kde(kernel, thresholds, data) {
  if (data.length == 0)
    return thresholds.map(t => [t,0]);
  return thresholds.map(t => [t, (d3.sum(data, d => kernel(t - d))/numberOfAllMutations)]); //get kde value at the bin boundaries
}


function epanechnikov(bw) {
  return function(v) {
    return Math.abs(v /= bw) <= 1 ? amplitude * (1 - v * v) / bw : 0;
  };
}


function filter_mutations(chr_number, data){ //filter dataset by mutations
  const data_chr = data.filter(x => x.chr === chr_number);
  let res = [];
  if (mutType == "point mutation"){
    for (let i = 0; i < data_chr.length; ++i) {
      res.push(parseInt(data_chr[i].pos));
    }
  }
  if (mutType == "genomic region"){
    for (let i = 0; i < data_chr.length; ++i) {
      res.push(parseInt(data_chr[i].start));
    }
  }
  return res
}


function _numberOfAllMutations(data){return(
data[0].length
)}



function _samples(data,getSamples)
{
  let res = []
  for (let i=0; i<data.length; i++){
    res = res.concat(["dataset"+String(i+1)])
    let tmp = getSamples(data[i])
    if (typeof tmp[0] == 'undefined')
      tmp = []
    res = res.concat(tmp)
  }
  return res
}


function getCheckboxOptions(color, samples){
  let res = [];
  for (let i=0; i<samples.length; i++){
    const x = {value: samples[i], label: samples[i], color: color[i] };
    res.push(x)
  }
  return res;
}


function getDensityOptions(samples){
  let res = ["all data"];
  for (let i=0; i<samples.length; i++){
    if (samples[i] != 'all data')
      res.push(samples[i]);
  }
  return res;
}



function checkbox(config = {}) {
  let {
    value: formValue, title, description, submit, options
  } = Array.isArray(config) ? {options: config} : config;
  options = options.map(
    o => (typeof o === "string" ? { value: o, label: o } : o)
  );
  const form = input({
    type: "checkbox",
    title,
    description,
    submit,
    getValue: input => {
      if (input.length)
        return Array.prototype.filter
          .call(input, i => i.checked)
          .map(i => i.value);
      return input.checked ? input.value : false;
    },
    form: html`
      <form>
        ${options.map(({ value, label, color}) => {
          const colorBlock = html`<span style="color: ${color}; font-weight: bold">&bull;</span>`
          const input = html`<input type=checkbox name=input ${
            (formValue || []).indexOf(value) > -1 ? "checked" : ""
          } style="vertical-align: baseline;" />`;
          input.setAttribute("value", value);
          const tag = html`<label style="display: inline-block; margin: 5px 10px 3px 0; font-size: 0.85em;">
           ${input}
           ${label}
           ${colorBlock}
          </label>`;
          return tag;
        })}
      </form>
    `
  });
  form.output.remove();
  return form;
}


function _d3(require){return(
require("d3@6")
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  function toString() { return this.url; }
  const fileAttachments = new Map([
    ["mouse_genome.txt", {url: new URL("./files/bba4899d367f46fab245abbede1fb61601ad628ae02f443acb589e7b26157d4182ae5eb104e1afc7a4caecddd95481db1e46eb4b7cb3b6e8a9495e7f98343f86.txt", import.meta.url), mimeType: "text/plain", toString}],
    ["human_genome.txt", {url: new URL("./files/5092fe2dcfc4d6066c3a27cfc2bcc8cbc6721d01a9786d65ed4e253399282baa22713fc18a8cc740c4bcce225fa71fa27f89ef8c033cc0dbf85e6ccb16890f53.txt", import.meta.url), mimeType: "text/plain", toString}]
  ]);
  main.builtin("FileAttachment", runtime.fileAttachments(name => fileAttachments.get(name)));
  main.variable(observer()).define(["md"], _1);
  main.variable(observer()).define(["md"], _2);
  main.variable(observer("viewof filesToUpload")).define("viewof filesToUpload", ["file"], _filesToUpload);
  main.variable(observer("filesToUpload")).define("filesToUpload", ["Generators", "viewof filesToUpload"], (G, _) => G.input(_));
  main.variable(observer("viewof genome")).define("viewof genome", ["select"], _genome);
  main.variable(observer("genome")).define("genome", ["Generators", "viewof genome"], (G, _) => G.input(_));
  main.variable(observer("viewof refGenome")).define("viewof refGenome", ["select"], _refGenome);
  main.variable(observer("refGenome")).define("refGenome", ["Generators", "viewof refGenome"], (G, _) => G.input(_));
  main.variable(observer("viewof mutType")).define("viewof mutType", ["select"], _mutType);
  main.variable(observer("mutType")).define("mutType", ["Generators", "viewof mutType"], (G, _) => G.input(_));
  main.variable(observer()).define(["md"], _7);
  main.variable(observer("viewof experiments")).define("viewof experiments", ["checkbox","getCheckboxOptions","color","samples"], _experiments);
  main.variable(observer("experiments")).define("experiments", ["Generators", "viewof experiments"], (G, _) => G.input(_));
  main.variable(observer()).define(["d3","genome","human_centromere","human_genome","mouse_centromere","mouse_genome","draw_chromosome","draw_mutations","data","color","experiments"], _9);
  main.variable(observer()).define(["md"], _10);
  main.variable(observer("viewof binSize")).define("viewof binSize", ["Inputs"], _binSize);
  main.variable(observer("binSize")).define("binSize", ["Generators", "viewof binSize"], (G, _) => G.input(_));
  main.variable(observer("viewof bandwidth")).define("viewof bandwidth", ["Inputs"], _bandwidth);
  main.variable(observer("bandwidth")).define("bandwidth", ["Generators", "viewof bandwidth"], (G, _) => G.input(_));
  main.variable(observer("viewof amplitude")).define("viewof amplitude", ["Inputs"], _amplitude);
  main.variable(observer("amplitude")).define("amplitude", ["Generators", "viewof amplitude"], (G, _) => G.input(_));
  main.variable(observer("viewof numMut")).define("viewof numMut", ["Inputs"], _numMut);
  main.variable(observer("numMut")).define("numMut", ["Generators", "viewof numMut"], (G, _) => G.input(_));
  main.variable(observer("viewof cluster")).define("viewof cluster", ["Inputs"], _cluster);
  main.variable(observer("cluster")).define("cluster", ["Generators", "viewof cluster"], (G, _) => G.input(_));
  main.variable(observer("viewof minClusterSize")).define("viewof minClusterSize", ["Inputs"], _minClusterSize);
  main.variable(observer("minClusterSize")).define("minClusterSize", ["Generators", "viewof minClusterSize"], (G, _) => G.input(_));
  main.variable(observer("viewof densData")).define("viewof densData", ["select","getDensityOptions","samples"], _densData);
  main.variable(observer("densData")).define("densData", ["Generators", "viewof densData"], (G, _) => G.input(_));
  main.variable(observer("viewof visions")).define("viewof visions", ["checkbox","numMut"], _visions);
  main.variable(observer("visions")).define("visions", ["Generators", "viewof visions"], (G, _) => G.input(_));
  main.variable(observer()).define(["all_diagramms","chrNumber"], _19);
  main.variable(observer()).define(["resultsCluster"], _20);
  main.variable(observer()).define(["generate_text_cluster","resultsCluster","md"], _21);
  main.variable(observer()).define(["md"], _22);
  main.variable(observer("viewof densData1")).define("viewof densData1", ["select","getDensityOptions","samples"], _densData1);
  main.variable(observer("densData1")).define("densData1", ["Generators", "viewof densData1"], (G, _) => G.input(_));
  main.variable(observer("viewof densData2")).define("viewof densData2", ["select","getDensityOptions","samples"], _densData2);
  main.variable(observer("densData2")).define("densData2", ["Generators", "viewof densData2"], (G, _) => G.input(_));
  main.variable(observer("comparison")).define("comparison", ["d3","genome","human_centromere","human_genome","mouse_centromere","mouse_genome","data","draw_diagram","getBinBoarders","densData1","densData2","draw_chromosome","draw_comparison_distribution"], _comparison);
  main.variable(observer()).define(["md"], _26);
  main.variable(observer()).define(["md"], _27);
  main.variable(observer("draw_chromosome")).define("draw_chromosome", _draw_chromosome);
  main.variable(observer("draw_mutations")).define("draw_mutations", ["draw_specific_mutations"], _draw_mutations);
  main.variable(observer("draw_specific_mutations")).define("draw_specific_mutations", ["mutType","length_mutation_centromere"], _draw_specific_mutations);
  main.variable(observer("length_mutation_centromere")).define("length_mutation_centromere", _length_mutation_centromere);
  main.variable(observer("draw_all_mutation_points")).define("draw_all_mutation_points", ["mutType"], _draw_all_mutation_points);
  main.variable(observer()).define(["md"], _33);
  main.define("initial chrNumber", _chrNumber);
  main.variable(observer("mutable chrNumber")).define("mutable chrNumber", ["Mutable", "initial chrNumber"], (M, _) => new M(_));
  main.variable(observer("chrNumber")).define("chrNumber", ["mutable chrNumber"], _ => _.generator);
  main.variable(observer("rgbValue")).define("rgbValue", _rgbValue);
  main.variable(observer("draw_comparison_distribution")).define("draw_comparison_distribution", ["getPolygoneOnChromosome","rgbValue"], _draw_comparison_distribution);
  main.variable(observer("draw_distribution")).define("draw_distribution", ["getPolygoneOnChromosome","d3","mutable chrNumber","generate_text"], _draw_distribution);
  main.variable(observer("getPolygoneOnChromosome")).define("getPolygoneOnChromosome", ["length_mutation_centromere"], _getPolygoneOnChromosome);
  main.variable(observer("define_cluster")).define("define_cluster", ["minClusterSize"], _define_cluster);
  main.variable(observer("draw_cluster")).define("draw_cluster", ["define_cluster","genome","linkGenomeBrowser"], _draw_cluster);
  main.variable(observer("mark_bins")).define("mark_bins", _mark_bins);
  main.variable(observer("linkGenomeBrowser")).define("linkGenomeBrowser", ["refGenome"], _linkGenomeBrowser);
  main.variable(observer("generate_text")).define("generate_text", _generate_text);
  main.variable(observer("generate_text_cluster")).define("generate_text_cluster", ["linkGenomeBrowser"], _generate_text_cluster);
  main.variable(observer("getBinBoarders")).define("getBinBoarders", _getBinBoarders);
  main.variable(observer("getMutationsGeneNumber")).define("getMutationsGeneNumber", _getMutationsGeneNumber);
  main.variable(observer("all_together")).define("all_together", ["d3","genome","human_centromere","human_genome","mouse_centromere","mouse_genome","getBinBoarders","all_diagramms","getMutationsGeneNumber","draw_chromosome","visions","draw_distribution","draw_all_mutation_points","mark_bins","numMut","draw_cluster","cluster"], _all_together);
  main.variable(observer("resultsCluster")).define("resultsCluster", ["densData","data","all_together","all_diagramms"], _resultsCluster);
  main.variable(observer()).define(["md"], _49);
  main.variable(observer("all_diagramms")).define("all_diagramms", ["genome","human_genome","mouse_genome","densData","data","draw_diagram"], _all_diagramms);
  main.variable(observer("draw_diagram")).define("draw_diagram", ["filter_mutations","d3","binSize","kde","epanechnikov","bandwidth","visions","cluster"], _draw_diagram);
  main.variable(observer("kde")).define("kde", ["d3","numberOfAllMutations"], _kde);
  main.variable(observer("epanechnikov")).define("epanechnikov", ["amplitude"], _epanechnikov);
  main.variable(observer("filter_mutations")).define("filter_mutations", ["mutType"], _filter_mutations);
  main.variable(observer()).define(["md"], _55);
  main.variable(observer("numberOfAllMutations")).define("numberOfAllMutations", ["data"], _numberOfAllMutations);
  main.variable(observer("data")).define("data", ["filesToUpload","convert_data","d3","Files","genome"], _data);
  main.variable(observer("getSamples")).define("getSamples", _getSamples);
  main.variable(observer("samples")).define("samples", ["data","getSamples"], _samples);
  main.variable(observer("getCheckboxOptions")).define("getCheckboxOptions", _getCheckboxOptions);
  main.variable(observer("getDensityOptions")).define("getDensityOptions", _getDensityOptions);
  main.variable(observer()).define(["md"], _62);
  main.variable(observer("convert_data")).define("convert_data", _convert_data);
  main.variable(observer("human_genome")).define("human_genome", ["d3","FileAttachment"], _human_genome);
  main.variable(observer("human_centromere")).define("human_centromere", ["human_genome"], _human_centromere);
  main.variable(observer("mouse_genome")).define("mouse_genome", ["d3","FileAttachment"], _mouse_genome);
  main.variable(observer("mouse_centromere")).define("mouse_centromere", ["mouse_genome"], _mouse_centromere);
  main.variable(observer()).define(["md"], _68);
  main.variable(observer("color")).define("color", _color);
  const child1 = runtime.module(define1);
  main.import("file", child1);
  const child2 = runtime.module(define1);
  main.import("input", child2);
  const child3 = runtime.module(define1);
  main.import("select", child3);
  main.variable(observer("checkbox")).define("checkbox", ["input","html"], _checkbox);
  main.variable(observer("d3")).define("d3", ["require"], _d3);
  main.variable(observer()).define(["md"], _75);
  main.variable(observer()).define(["md"], _76);
  main.variable(observer()).define(["md"], _77);
  return main;
}
