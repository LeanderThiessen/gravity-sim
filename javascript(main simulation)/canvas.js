const canvas = document.querySelector('canvas')
const c = canvas.getContext('2d')
c.font = "30px Arial";
canvas.width = innerWidth
canvas.height = innerHeight

const mouse = {
  x: innerWidth / 2,
  y: innerHeight / 2
}

// Event Listeners
addEventListener('mousemove', (event) => {
  mouse.x = event.clientX
  mouse.y = event.clientY
})

addEventListener('resize', () => {
  canvas.width = innerWidth
  canvas.height = innerHeight

  init()
})

// Objects
class Object {
  constructor(x, y, dx, dy, radius, color) {
    this.x = x
    this.y = y
    this.x_Pos = []
    this.y_Pos = []
    this.dx = dx
    this.dy = dy
    this.radius = radius
    this.color = color
    this.verletList = []
    this.coll_record = this
  }

  update() {
   
    this.x_Pos.push(this.x)
    this.y_Pos.push(this.y)
    
    let F_x = 0
    let F_y = 0

    //Sum over all the forces acting on 'this' from objects in Verlet List
    this.verletList.forEach(object => {
      let x_diff = this.x - object.x
      let y_diff = this.y - object.y
      
      if (boundary_cond === 'periodic'){
        //extend interaction over pbc
        if (x_diff > r_cutoff){
          x_diff -= canvas.width
        }
        if (x_diff < -r_cutoff){
          x_diff += canvas.width
        }

        if (y_diff > r_cutoff){
          y_diff -= canvas.height
        }
        if (y_diff < -r_cutoff){
          y_diff += canvas.height
        }
      }

      let r = Math.sqrt(x_diff**2+y_diff**2)
      //F shouldn't be calculated for smaller distances than R, velocities would blow up
      if (r > G_cutoff){
        F_x += -(G/r**3/*-K/r**p*/)*x_diff
        F_y += -(G/r**3/*-K/r**p*/)*y_diff
      }
    })


    //Update positions
    let new_x = this.x + this.dx*dt + F_x*dt**2
    let new_y = this.y + this.dy*dt + F_y*dt**2

    let F_new_x = 0
    let F_new_y = 0
    //compute forces on 'this' at updated positions
    this.verletList.forEach(object => {
      let x_diff = new_x - object.x
      let y_diff = new_y - object.y

      if (boundary_cond === 'periodic'){
        //extend interaction over pbc
        if (x_diff > r_cutoff){
          x_diff -= canvas.width
        }
        if (x_diff < -r_cutoff){
          x_diff += canvas.width
        }

        if (y_diff > r_cutoff){
          y_diff -= canvas.height
        }
        if (y_diff < -r_cutoff){
          y_diff += canvas.height
        }
      }
      let r = Math.sqrt(x_diff**2+y_diff**2)
      
      if (r > G_cutoff){
        F_new_x += -(G/r**3/*-K/r**p*/)*x_diff 
        F_new_y += -(G/r**3/*-K/r**p*/)*y_diff
      }
    })


    let new_v_x = this.dx + (F_x + F_new_x)/2*dt
    let new_v_y = this.dy + (F_y + F_new_y)/2*dt + g*dt

    this.x = new_x
    this.y = new_y

    this.dx = new_v_x
    this.dy = new_v_y

    
    if (boundary_cond === 'periodic'){
      if (new_x < 0){
        this.x = canvas.width
      }
      if (new_x > canvas.width){
        this.x = 0//new_x % canvas.width
      }

      if (new_y < 0){
        this.y = canvas.height
      }
      if (new_y > canvas.height){
        this.y = 0//new_y % canvas.height
      }
    }

    
    else if (boundary_cond === 'fixed'){
      let delta = 10
      if (new_x < delta){
        this.x = delta
        this.dx *= -1
        this.dx *= fric
        this.dy *= fric
      }
      if (new_x > canvas.width-delta){
        this.x = canvas.width-delta
        this.dx *= -1
        this.dx *= fric
        this.dy *= fric
      }

      if (new_y < delta){
        this.y = delta
        this.dy *= -1
        this.dx *= fric
        this.dy *= fric
      }
      if (new_y > canvas.height-delta){
        this.y = canvas.height-delta
        this.dy *= -1
        this.dx *= fric
        this.dy *= fric
      }
    }
  }

  updateVerletList() {
    this.verletList = []
    for (let i = 0; i < N; i++){
      let push = false
      if (objects[i] != this){

        let dist_squared = (objects[i].x-this.x)**2+(objects[i].y-this.y)**2 
        
        if (dist_squared <= r_cutoff**2 ){
          push = true
        }

        else  {
          //x-direction periodic boundary conditions (right->left)
          let x_new = (objects[i].x - canvas.width)
          let y_new = objects[i].y 
          dist_squared = (x_new-this.x)**2+(y_new-this.y)**2 

          if (dist_squared <= r_cutoff**2 ){
            push = true
          }

          //x-direction periodic boundary conditions (left->right)
          x_new = (objects[i].x + canvas.width)
          y_new = objects[i].y 
          dist_squared = (x_new-this.x)**2+(y_new-this.y)**2 

          if (dist_squared <= r_cutoff**2 ){
            push = true
          }

          //y-direction periodic boundary conditions (up->down)
          x_new = objects[i].x 
          y_new = objects[i].y - canvas.height
          dist_squared = (x_new-this.x)**2+(y_new-this.y)**2 

          if (dist_squared <= r_cutoff**2 ){
            push = true
          }

          //y-direction periodic boundary conditions (down->up)
          x_new = objects[i].x 
          y_new = objects[i].y + canvas.height
          dist_squared = (x_new-this.x)**2+(y_new-this.y)**2 

          if (dist_squared <= r_cutoff**2 ){
            push = true
          }
        }
      }
      if (push === true){
        this.verletList.push(objects[i])
      }
    }
  }
}

class simTime {
  constructor (){
    this.t = 0
  }
  tick(){
    this.t ++;
  }
}

function collisionCheck(objects){
  //loop through all unique pairs of different objects
  for (let i = 0; i < N; i++){
    let currentVerletList = objects[i].verletList
    
    for (let j = 0; j < currentVerletList.length; j++){

      //only compute the collision if object i hasn't yet collided with object j in the same update loop/frame
      if (objects[i].coll_record != currentVerletList[j]){
        //cumpute distance vector between balls
        let d_x = objects[i].x-currentVerletList[j].x
        let d_y = objects[i].y-currentVerletList[j].y
          
        let d_abs_squared = d_x**2+d_y**2
          
        //check whether balls are overlapping
        if (d_abs_squared < (2*objects[i].radius)**2){
            
          let x_vel = currentVerletList[j].dx - objects[i].dx
          let y_vel = currentVerletList[j].dy - objects[i].dy

          let vel_dot = d_x*x_vel + d_y*y_vel
            
          //only perform collision when balls are moving towards each other (otherwise they get stuck)
          if (vel_dot > 0){
            let d_abs = Math.sqrt(d_abs_squared)
            d_x = d_x/d_abs
            d_y = d_y/d_abs
            let a = d_x*(objects[i].dx - currentVerletList[j].dx) + d_y*(objects[i].dy - currentVerletList[j].dy);
              
            //update new velocities after collision
            objects[i].dx -= a*d_x
            objects[i].dy -= a*d_y           


            currentVerletList[j].dx += a*d_x
            currentVerletList[j].dy += a*d_y

            //everytime the balls overlap, reset positions such that they don't overlap anymore
            currentVerletList[j].x = objects[i].x - 2*objects[i].radius*d_x
            currentVerletList[j].y = objects[i].y - 2*objects[i].radius*d_y

            currentVerletList[j].coll_record = objects[i]
          }
        }
      }
    }
  }
}


function init() {
  if (initial === 'Random'){
    let j = 0;

    for (let i = 0; i<N; i++){
      let v_x = 0*(2*Math.random()-1)
      let v_y = 0*(2*Math.random()-1)

      let delta = 1//7*radius
      let N_W = N/canvas.height*delta
      let N_H = N/canvas.width*delta

      let x_0 = canvas.width*Math.random()//((i+1)*delta) % (canvas.width-radius)
      let y_0 = canvas.height*Math.random()//delta * Math.floor(((i+1)*delta) / (canvas.width-radius))
      
      objects.push(new Object(x_0,y_0,v_x,v_y, radius,'blue'))
    }
  }

  if (initial === 'two_circles'){
    let R_max = 30
    let lam = 15
    for (let i = 0; i<N/2; i++){
      let phi = Math.random()*2*Math.PI
      let R = R_max*Math.random()

      let v_x = -lam*R*Math.sin(phi) + 800
      let v_y = lam*R*Math.cos(phi)

      let x_0 = R*Math.cos(phi)+0.4*canvas.width
      let y_0 = R*Math.sin(phi)+0.5*canvas.height
      
      objects.push(new Object(x_0,y_0,v_x,v_y, radius,'red'))
    }

    R_max = 30
    
    for (let i = 0; i<N; i++){
      let phi = Math.random()*2*Math.PI
      let R = R_max*Math.random()

      let v_x = -lam*R*Math.sin(phi) -800
      let v_y = lam*R*Math.cos(phi)

      let x_0 = R*Math.cos(phi)+0.6*canvas.width
      let y_0 = R*Math.sin(phi)+0.5*canvas.height
      



      objects.push(new Object(x_0,y_0,v_x,v_y, radius,'red'))
    }

  }

  if (initial === 'one_circle'){
    /*
    const simLength = 200 //number of computed timesteps
    const dt = 0.0018 //timestep
    const N = 4000 //number of particles
    const radius = 0.8 //radius of particles (for collisions)
    const G = 6000  //gravitational strength
    const G_cutoff = 2*radius //cutoff distance of gravitational force (for close contact), ignore forces below this scale
    const r_scale = 1 //scales the drawn radius of particles by r_scale
    */
    let R_max = 100
    
    for (let i = 0; i<N; i++){
      let phi = Math.random()*2*Math.PI
      let R = R_max*Math.random()

      //let lam = 38*(1-R/R_max)**2 + 5.2
      //let v_x = -R*Math.sin(phi)*lam
      //let v_y = R*Math.cos(phi)*lam

      let v_x = 0
      let v_y = 0

      let x_0 = R*Math.cos(phi)+0.5*canvas.width
      let y_0 = R*Math.sin(phi)+0.75*canvas.height
      

      objects.push(new Object(x_0,y_0,v_x,v_y, radius,'red'))
    }
  }
}

function integrateEom(){
  for (let t = 0; t < simLength; t++){
    //Collision-Check
    if (withCollision === true){
      collisionCheck(objects)
    }

    //Update Verlet List
    objects.forEach(object => {
      object.coll_record = this
      if ((t % n_list) === 0){
        object.updateVerletList()
      }
    })

    //Position and velocity update (Velocity Verlet)
    objects.forEach(object => {
      object.update()
    })
  }
}

//Animation
function animate() {
  requestAnimationFrame(animate)

  if (withTrace == false){
    c.clearRect(0, 0, canvas.width, canvas.height)
  }

  objects.forEach(object => {
    let x = object.x_Pos[elapsedTime.t % simLength]
    let y = object.y_Pos[elapsedTime.t % simLength]

    c.beginPath()
    c.arc(x, y, r_scale*object.radius, 0, Math.PI * 2, false)
    c.fillStyle = object.color
    c.fill()
    c.closePath()
  })
  
  elapsedTime.tick()
  
  if (elapsedTime.t % simLength === 0){
    c.clearRect(0, 0, canvas.width, canvas.height)
  } 
}


//Variables
let objects = [];
let elapsedTime = new simTime;
const simLength = 500 //number of computed timesteps
const dt = 0.001 //timestep
const N = 1200 //number of particles
const radius = 2//radius of particles (for collisions)
const G = 500000  //gravitational strength
const fric = 0.999
const G_cutoff = 2*radius //cutoff distance of gravitational force (for close contact), ignore forces below this scale
const r_scale = 1 //scales the drawn radius of particles by r_scale
const g = 0 //gravitational pull to the bottom
const n_list = 10 //number of timesteps between Verlet list update
const colors = ['blue'] 
const r_cutoff = canvas.width/4//30*radius//canvas.width/2 //range cutoff for gravitational force/ Verlet list
const initial = 'Random' //initial configuration of particle ('Random','Circle')
const boundary_cond = 'periodic' //'periodic','fixed'
const withCollision = true //turn collisions off/on
const withTrace = false //if true, window is not cleared between animation timesteps


//Run
console.log("Initialize")
init()
console.log("Integrate eom...")
integrateEom()
console.log("Prepare animation...")
animate()

