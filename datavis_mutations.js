
function draw_svg()
{
  /*important constants*/
  const width = $("#datavis").width(), //width of the image 
    height = width*1.2, //height of the image
    c_width = width/50, //width of one chromosome
    distance_c = 2*c_width,  //distance between two chromosomes 
    c_factor = 320000, //scaling factor for the chromosomes
    centromere_x = c_width/5, //distance for centromere in x-direction
    centromere_y = c_width/2; //distance where centromere starts/ends in y-direction 

  const zoom = d3.zoom()
      .scaleExtent([1, 20])
      .on("zoom", zoomed),
    svg = d3
      .select("#datavis")
      .attr("viewBox", [0, 0,width, height]),
    g = svg.append("g");

  let centromere, data_genome;
  if (sessionStorage.getItem("genomeType") == "human"){
    centromere = human_centromere;
    data_genome = human_genome;
  } else if (sessionStorage.getItem("genomeType") == "mouse"){
    centromere = mouse_centromere;
    data_genome = mouse_genome;
  }

  /*draw chromosomes depending on genome type*/
  let centromere_d = draw_chromosome(g, data_genome, centromere, distance_c, c_width, c_factor, centromere_x, centromere_y);

  /*get all checked Checkboxes as an array*/
  var allChecked = getCheckedCheckboxes( $('.filled-in').get());
  
  /*draw mutations of dataset 1*/
  var data1 = JSON.parse(sessionStorage.getItem("data1"));
  if(data1.length != 0)
    draw_mutations(data1, color1, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y, allChecked, data_genome, sessionStorage.getItem("datatype1"), 1, sessionStorage.getItem("value"), parseFloat(sessionStorage.getItem("startValue")), parseFloat(sessionStorage.getItem("endValue")));
  
  /*draw mutations of dataset2*/
  var data2 = JSON.parse(sessionStorage.getItem("data2"));
  if(data2.length != 0)
    draw_mutations(data2, color2, g, distance_c, c_width, c_factor, centromere_d, centromere_x, centromere_y, allChecked, data_genome, sessionStorage.getItem("datatype2"), 2, sessionStorage.getItem("value"), parseFloat(sessionStorage.getItem("startValue")), parseFloat(sessionStorage.getItem("endValue")));

  /*allows zooming*/
  svg.call(zoom);
  
  function zoomed() {
    const {transform} = d3.event;
    g.attr("transform", transform);
    g.attr("stroke-width", 1 / transform.k);
  }
  
  return svg.node();
}

function getCheckedCheckboxes(checkboxes){
  var checked = [];
  for(var i=0; i<checkboxes.length; i++){
    if (checkboxes[i].checked){
      checked.push(checkboxes[i].id)
    }
  }
  return checked;
}


function draw_chromosome(g, data_genome, centromere, distance_c, c_width, c_factor, centromere_x, centromere_y){
  /*
    draw chromosomes of given genome
    returns the positions of the drwan centromere
  */
  const centromere_d = []; //array stores positions of all centromere  
  /*fill array with positions of all centromere*/
  for (let i = 0; i < data_genome.length; ++i) {
    //use average of centromere in the data for visualised centromere
    centromere_d.push((parseInt(centromere[i][0]) + parseInt(centromere[i][1]))/(2*c_factor)) 
  }

  /*Draw Chromosome*/
  for (let i = 0; i < data_genome.length; ++i) {
    const c = centromere_d[i]; //coordinate of centromere
    let beginnCentromere = 0;
    //start of centromer has to fit inside chromosome boundaries
    if (c-centromere_y < 0)
      beginnCentromere = 0;
    else
      beginnCentromere = c-centromere_y
    const chromosome = g.append("polyline")
    .attr("points", String(i*distance_c) + " " + String(0) + ", " //top of chromosome (left side)
          + String(i*distance_c) + " " + String(beginnCentromere) + ", " //start of centromere
          + String(i*distance_c + centromere_x) + " " + String(c) + ", " //centromere
          + String(i*distance_c) + " " + String(c + centromere_y) + ", " //end of centromere
          + String(i*distance_c) + " " + String(data_genome[i].size/c_factor) + ", " //bottom of chromosome (left side)
          + String(i*distance_c + c_width) + " " + String(data_genome[i].size/c_factor) + ", " //bottom (right side)
          + String(i*distance_c + c_width) + " " + String(c + centromere_y) + ", " //end centromere
          + String(i*distance_c + c_width - centromere_x) + " " + String(c) + ", " //centromere
          + String(i*distance_c + c_width) + " " + String(beginnCentromere) + ", " //start centromere
          + String(i*distance_c + c_width) + " " + String(0) + ", " //top (right side)
          + String(i*distance_c) + " " + String(0)) //top (left side) for closed line
    .style("fill", "#fff")
    .style("stroke", "black");
    //show ID when moving mouse over chromosome
    chromosome.append("title") 
    .text("chr " + data_genome[i].ID);
    //label all chromosome with ID
    g.append("text")
      .attr("x", i*distance_c + c_width/2)             
      .attr("y", 15 + data_genome[i].size/c_factor)
      .attr("text-anchor", "middle")  
      .style("font-size", "12px") 
      .text(String(data_genome[i].ID));
  }
  return centromere_d;
}


function draw_mutations(data, colors, svg, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y, layers, genome, mutType, dataNum, value, start, end){ 
  /*
    visualise mutations from selected experiments
  */

  var sample1 = JSON.parse(sessionStorage.getItem('samples1'));
  var sample2 = JSON.parse(sessionStorage.getItem('samples2'));

  for (let i = 0; i < layers.length; ++i) {
    if (layers[i] == "Dataset1"){
      if (dataNum == 1)
        draw_specific_mutations(data, svg, colors[0], distance_c, c_width,c_factor, centromere, centromere_x, centromere_y, genome, mutType, value, start, end);
    }else if (layers[i] == "Dataset2"){
      if (dataNum == 2)
        draw_specific_mutations(data, svg, colors[0], distance_c, c_width,c_factor, centromere, centromere_x, centromere_y, genome, mutType, value, start, end);
    }else{
      //filter data by experiment
      let ind1 = sample1.indexOf(layers[i]);
      let ind2 = sample2.indexOf(layers[i]);

      if (ind1 != -1){
        let experiment = []
        for (let j=0; j<data.length; j++){
          experiment = data.filter(x => x.sample === layers[i]); 
        }
        //draw mutations of experiment
        draw_specific_mutations(experiment, svg, colors[ind1] , distance_c, c_width,c_factor, centromere, centromere_x, centromere_y, genome, mutType, value, start, end); 
      } else if (ind2 != -1){
        let experiment = []
        for (let j=0; j<data.length; j++)
          experiment = data.filter(x => x.sample === layers[i]); 
        //draw mutations of experiment
        draw_specific_mutations(experiment, svg, colors[ind2] , distance_c, c_width,c_factor, centromere, centromere_x, centromere_y, genome, mutType, value, start, end); 
      }
    }
  }
}


function draw_specific_mutations(experement_data, g, color, distance_c, c_width, c_factor, centromere, centromere_x, centromere_y, genome, mutType, value, start, end){
  /*
    visualise mutations on the chromosomes for one experiment
  */
  for (let i = 0; i < experement_data.length; ++i) { 
    const chr = parseInt(experement_data[i].chr); //chromosome of the mutation

    //change color dependent on value
    var col;
    if(value == "true"){
      var myColor = d3.scaleLinear().domain([start,end]).range([color, "black"]); //geht nur für explizite Namen...
      col = myColor(experement_data[i].value);
      console.log(col)
    }
    else{
      col = color;
    }
    
    //To-Do: mousemove over mutation --> get value

    if((chr-1) < genome.length){ //only mutations on existing chromosomes is visualised
      //data outside the chromosomes is clipped away --> error message/notification
      if (mutType == "pm"){ //draw point mutation as a line
        let position = experement_data[i].pos/c_factor;
        if (position < genome[chr-1].size/c_factor){ //only mutations on the length of the chromosome is visualised
          if ((position > (centromere[chr-1] - centromere_y))  && (position < (centromere[chr-1] + centromere_y))){   
            //mutation in the area of the centromere
            //change length of the mutation line
            let distance_m = length_mutation_centromere(centromere[chr-1], centromere_x, centromere_y, position) 
            //draw mutation as a line
            let mut = g.append("line")
              .attr("x1", (chr-1)*distance_c + distance_m)
              .attr("x2", (chr-1)*distance_c + c_width - distance_m)
              .attr("y1", position) 
              .attr("y2", position) 
              .attr("stroke", col);
            mut.append("title").text(experement_data[i].value);
          } else { 
            //mutation not in the area of the centromere
            //draw mutation as a line with standard length
            let mut = g.append("line")
              .attr("x1", (chr-1)*distance_c)
              .attr("x2", (chr-1)*distance_c + c_width)
              .attr("y1", position) 
              .attr("y2", position) 
              .attr("stroke", col);
            mut.append("title").text(experement_data[i].value);
          }
        }
      } else if (mutType == "gr") { //draw genomic region as a polygon
        let positionStart = experement_data[i].start/c_factor;
        let positionEnd = experement_data[i].end/c_factor;
        if (positionEnd < genome[chr-1].size/c_factor){ //only genomic regions on the chromosomes are visualised
          //müssen noch am Centromer angepasst werden!
          let distance_left_start = (chr-1)*distance_c;
          let distance_right_start = (chr-1)*distance_c;
          let distance_left_end = (chr-1)*distance_c + c_width;
          let distance_right_end = (chr-1)*distance_c + c_width;
          if ((positionStart > (centromere[chr-1] - centromere_y))  && (positionStart < (centromere[chr-1] + centromere_y))){
            let tmp = length_mutation_centromere(centromere[chr-1], centromere_x, centromere_y, positionStart)
            distance_left_start += tmp;
            distance_right_start -= tmp;
          }
          if ((positionEnd > (centromere[chr-1] - centromere_y))  && (positionEnd < (centromere[chr-1] + centromere_y))){
            let tmp = length_mutation_centromere(centromere[chr-1], centromere_x, centromere_y, positionEnd)
            distance_left_end += tmp;
            distance_right_end -= tmp;
          }
          let polygon = String(distance_left_start) + " " + String(positionStart) + ", " + 
            String(distance_left_end) + " " + String(positionEnd) + ", " +
            String(distance_right_end) + " " + String(positionEnd) + ", " + 
            String(distance_right_start) + " " + String(positionStart)+ ", " + 
            String(distance_left_start) + " " + String(positionStart);
          let mut = g.append("polyline")
            .attr("points", polygon) 
            .style("stroke", col)
            .style("fill", col)
          mut.append("title").text(experement_data[i].value);
        }
      }
    }
  }
}


function length_mutation_centromere(centromere, centromere_x, centromere_y, mutation_pos){ 
  /*
    change the length of the mutation line at the area of the centromere
    return the change of the length
  */
  let slope = centromere_y/centromere_x,
    centromere_bottom = centromere + centromere_y,
    centromere_top = centromere - centromere_y;
  if (mutation_pos < centromere) {
    let diff = centromere - mutation_pos
    return (centromere_bottom - (centromere + diff))/slope;
  } else {
    return (centromere_bottom - mutation_pos)/slope;
  }
}


function draw_all_mutation_points(data, g, distance_c, c_width, c_factor, genome){
  /*
    draw all mutations as half red lines, used to show it on top of the density
  	also genomic regions are drawn as lines
  */ 
  for (let i = 0; i < data.length; ++i) {
    var position = position = data[i].pos/c_factor;
    if(isNaN(position))
      position = data[i].start/c_factor;
    let chr = parseInt(data[i].chr);
    if(chr <= genome.length){
      if (position < genome[chr-1].size/c_factor){
        g.append("line")
          .attr("x1", (chr-1)*distance_c + 0.5 * c_width)
          .attr("x2", (chr-1)*distance_c + c_width)
          .attr("y1", position) 
          .attr("y2", position) 
          .attr("stroke", "red");
      }
    }
  }
}



