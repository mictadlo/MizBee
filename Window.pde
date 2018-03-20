float TRIANGLE_HEIGHT = 25.0;
float HALF_TRIANGLE_HEIGHT = 0.5*TRIANGLE_HEIGHT;
float TRIANGLE_VERTEX_AXIS_OFFSET = TRIANGLE_HEIGHT/sqrt(3.0);

class Window 
{
  float _deg_start, _deg_end, _deg_length;
  float _circle_r, _chrom_w;

  float[] _shutter_start, _shutter_end; // centers of the shutter trianlges

  Window( float ds, float de, float r, float w )
  {
    _deg_start = ds;
    _deg_end = de;
    _deg_length = de - ds;

    _shutter_start = new float[2];
    _shutter_end = new float[2];

    setRadiusAndWidth( r, w );
  }

  void setRadiusAndWidth( float r, float w ) 
  { 
    _circle_r = r; 
    _chrom_w = w; 

    // update the shutter positions 
    setShutterPositions();   
  }

  float startDegree() 
  { 
    return _deg_start; 
  }

  float endDegree()
  { 
    return _deg_end; 
  }

  void setStartAndEndDegrees( float s, float e ) {
    _deg_start = s;
    _deg_end = e;
    setShutterPositions();
  }

  void setShutterPositions() {
    _shutter_start[0] = (_circle_r+_chrom_w+HALF_TRIANGLE_HEIGHT)*cos(_deg_start);
    _shutter_start[1] = (_circle_r+_chrom_w+HALF_TRIANGLE_HEIGHT)*sin(_deg_start);

    _shutter_end[0] = (_circle_r+_chrom_w+HALF_TRIANGLE_HEIGHT)*cos(_deg_end);
    _shutter_end[1] = (_circle_r+_chrom_w+HALF_TRIANGLE_HEIGHT)*sin(_deg_end);
  }


  boolean overStartShutter( float x, float y ) {
    float x_diff = abs((float)x - _shutter_start[0]);
    float y_diff = abs((float)y - _shutter_start[1]);

    if ( (x_diff <= HALF_TRIANGLE_HEIGHT) &&
      (y_diff <= HALF_TRIANGLE_HEIGHT) ) {
      return true; 
    }
    else {      
      return false;
    }
  }

  boolean overEndShutter( float x, float y ) {
    float x_diff = abs((float)x - _shutter_end[0]);
    float y_diff = abs((float)y - _shutter_end[1]);

    if ( (x_diff <= HALF_TRIANGLE_HEIGHT) &&
      (y_diff <= HALF_TRIANGLE_HEIGHT) ) {
      return true; 
    }
    else {      
      return false;
    }
  }

  void updateShutter( float x, float y, boolean start_shutter, boolean end_shutter, float start_deg, float stop_deg ) {
    // get the degree of the new position
    float deg = (float)Math.atan2( (double)y,(double)(x+1.0e-6) );     

    if ( start_shutter ) {
      _deg_start = min(max(deg, start_deg), _deg_end);
    }
    else {
      _deg_end = max(min(deg, stop_deg), _deg_start);
    }

    setShutterPositions();
  }

  void render() {

    smooth();
    noStroke();

    // draw the window on the chromosome    
    fill( 255 );     

    if ( _deg_start <= _deg_end ) {
      beginShape( QUAD_STRIP );
      {
        for ( float a = _deg_start; a < _deg_end; a += DEG_STEPSIZE )
        {
          vertex( _circle_r*cos(a), _circle_r*sin(a) );
          vertex( (_circle_r+_chrom_w)*cos(a), (_circle_r+_chrom_w)*sin(a) );
        }

        vertex( _circle_r*cos(_deg_end), _circle_r*sin(_deg_end) );
        vertex( (_circle_r+_chrom_w)*cos(_deg_end), (_circle_r+_chrom_w)*sin(_deg_end) );
      }
      endShape();
    }
    else {
      beginShape( QUAD_STRIP );
      {
        for ( float a = _deg_start; a < PI; a += DEG_STEPSIZE )
        {
          vertex( _circle_r*cos(a), _circle_r*sin(a) );
          vertex( (_circle_r+_chrom_w)*cos(a), (_circle_r+_chrom_w)*sin(a) );
        }
        for ( float a = -PI; a < _deg_end; a += DEG_STEPSIZE )
        {
          vertex( _circle_r*cos(a), _circle_r*sin(a) );
          vertex( (_circle_r+_chrom_w)*cos(a), (_circle_r+_chrom_w)*sin(a) );
        }

        vertex( _circle_r*cos(_deg_end), _circle_r*sin(_deg_end) );
        vertex( (_circle_r+_chrom_w)*cos(_deg_end), (_circle_r+_chrom_w)*sin(_deg_end) );
      }
      endShape();
    }

    // draw the triangle markers
    fill( color(0) );
    strokeWeight(1.0);
    stroke( _chrom_color );

    pushMatrix();
    rotate( _deg_start );
    triangle( _circle_r+_chrom_w, 0.0, 
    _circle_r+_chrom_w+TRIANGLE_HEIGHT, -TRIANGLE_VERTEX_AXIS_OFFSET,
    _circle_r+_chrom_w+TRIANGLE_HEIGHT, TRIANGLE_VERTEX_AXIS_OFFSET );  
    popMatrix();   

    pushMatrix();
    rotate( _deg_end );
    triangle( _circle_r+_chrom_w, 0.0, 
    _circle_r+_chrom_w+TRIANGLE_HEIGHT, -TRIANGLE_VERTEX_AXIS_OFFSET,
    _circle_r+_chrom_w+TRIANGLE_HEIGHT, TRIANGLE_VERTEX_AXIS_OFFSET );  
    popMatrix();        
  }

}

