use plotters::prelude::*;
use crate::R;

#[allow(non_upper_case_globals)]
const Lx: u32 = 640;
#[allow(non_upper_case_globals)]
const Ly: u32 = 480;

#[allow(clippy::too_many_arguments)]
pub fn plot_1d(x: &[Vec<R>], y: &[Vec<R>], fname: &str,
        title: &str, xlabel: &str, ylabel: &str, labels: &[String],
        col1: (u8, u8, u8), col2: (u8, u8, u8)) 
    -> Result<(), Box<dyn std::error::Error>> 
{
     
    // check that we have the right number of labels
    if labels.len() != y.len() {
        return Err(TextError::new_box(format!("Different numbers of labels and lines ({} and {})",
                            labels.len(), y.len())));
    }
 
    let fname = format!("{}.svg", fname);
    let root = SVGBackend::new(&fname, (Lx, Ly)).into_drawing_area();
        
    // find the minimum and maximum of x and y
    let (x_min, x_max) = find_min_max(x)?;
    let (y_min, y_max) = find_min_max(y)?;
     
    root.fill(&WHITE)?;
    let root = root.margin(20, 20, 20, 20);
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 20).into_font())
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

    chart
        .configure_mesh()
        .x_labels(10)
        .y_labels(10)
        .x_desc(xlabel)
        .y_desc(ylabel)
        .axis_desc_style(("sans-serif", 15))
        .draw()?;
    
    let n_plots = y.len();
    for j in 0..n_plots {

        let x = x[j % x.len()].clone();
        let y = y[j].clone();

        // check that x and y have the same size
        if x.len() != y.len() {
            return Err(TextError::new_box(format!("The x and y vectors have different sizes ({} and {})",
                                x.len(), y.len())));
        }
        
        let ep = if n_plots > 1 {
            (j as f64) / (n_plots - 1) as f64
        } else {
            0.
        };
        let color = RGBColor(
            ((col2.0 as f64) * ep + (col1.0 as f64) * (1.-ep)) as u8, 
            ((col2.1 as f64) * ep + (col1.1 as f64) * (1.-ep)) as u8, 
            ((col2.2 as f64) * ep + (col1.2 as f64) * (1.-ep)) as u8,
        );
         
        // build the vector of coordinates
        let mut z = Vec::<(R,R)>::new();
        for i in 0..x.len() {
            z.push((x[i], y[i]));
        }
        let z = z;

        chart.draw_series(LineSeries::new(z.clone(), &color))?
            .label(&labels[j])
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &color));
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.5))
        .border_style(&BLACK.mix(0.5))
        .position(plotters::chart::SeriesLabelPosition::LowerRight)
        .draw()?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn plot_2d(z: &[R], fname: &str, 
        x_min: R, x_max: R, y_min: R, y_max: R,
        nx: usize, ny: usize, 
        show_axes: bool)
    -> Result<(), Box<dyn std::error::Error>> 
{
     
    // check that z has length nx * ny
    if z.len() != nx * ny {
        return Err(TextError::new_box(format!("Expected z.len() = nx * ny, but got {} and {}",
                                              z.len(), nx * ny)));
    }
 
    let fname = format!("{}.png", fname);
    let (lx, ly) = (Ly, ((Ly as R) * (y_max-y_min) / (x_max-x_min)) as u32);
    let root = BitMapBackend::new(&fname, (lx, ly)).into_drawing_area();
        
    // find the minimum and maximum of x, y, and z
    let (z_min, z_max) = find_min_max(&[z.to_vec()])?;
     
    root.fill(&WHITE)?;
    let root = if show_axes { 
        root.margin(20, 20, 20, 20)
    } else {
        root.margin(0, 0, 0, 0)
    };
    let mut chart = if show_axes { 
        ChartBuilder::on(&root)
            .x_label_area_size(20)
            .y_label_area_size(20)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)?
    } else {
        ChartBuilder::on(&root)
            .x_label_area_size(0)
            .y_label_area_size(0)
            .build_cartesian_2d(x_min..x_max, y_min..y_max)?
    };

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;
   
    let dx = (x_max-x_min) / (nx as R - 1.);
    let dy = (y_max-y_min) / (ny as R - 1.);
    let plotting_area = chart.plotting_area();
    let range = plotting_area.get_pixel_range();
    let (pw, ph) = (range.0.end - range.0.start, range.1.end - range.1.start);
    let (xr, yr) = (chart.x_range(), chart.y_range());
    for a in 0..=ph {
        for b in 0..=pw {
            let xr = xr.start + (b as R / (pw as R)) * (xr.end - xr.start);
            let yr = yr.start + (a as R / (ph as R)) * (yr.end - yr.start);
            let ir = (xr-x_min) / dx;
            let jr = (yr-y_min) / dy;
            let i = ir.floor();
            let j = jr.floor();
            let ep_x = ir - i;
            let ep_y = jr - j;
            let i = i as usize;
            let j = j as usize;
            let zr = if (i+1 < nx) && (j+1 < ny) {
                (1.-ep_x) * (1.-ep_y) * z[i*ny+j] 
                + (1.-ep_x) * ep_y * z[i*ny+j+1]
                + ep_x * (1. - ep_y) * z[(i+1)*ny+j]
                + ep_x * ep_y * z[(i+1)*ny+j+1]
            } else if i+1 < nx {
                (1.-ep_x) * z[i*ny+j] 
                + ep_x * z[(i+1)*ny+j]
            } else if j+1 < ny {
                (1.-ep_y) * z[i*ny+j] 
                + ep_y * z[i*ny+j+1]
            } else {
                z[i*ny+j]
            };
            let c = (zr-z_min)/(z_max-z_min);
            if c < 0. {
                plotting_area.draw_pixel((xr, yr), &HSLColor(0., 1.0, 0.5))?;
            } else if c > 1. {
                plotting_area.draw_pixel((xr, yr), &HSLColor(1., 1.0, 0.5))?;
            } else {
                plotting_area.draw_pixel((xr, yr), &HSLColor(c, 1.0, 0.5))?;
            }
        }
    }

    #[allow(clippy::expect_fun_call)]
    root.present().expect(&format!("Unable to write fo file {}", fname));

    Ok(())
}


fn find_min_max<T: PartialOrd + Copy>(x: &[Vec<T>]) -> Result<(T,T), TextError> {
    
    let mut x_min: Option<T> = None;
    let mut x_max: Option<T> = None;
    
    for x in x.iter() {
        for &x in x.iter() {
            if let (Some(x_min_c), Some(x_max_c)) = (x_min, x_max) {
               if x < x_min_c { x_min = Some(x) };
                if x > x_max_c { x_max = Some(x) };
            } else {
                x_min = Some(x);
                x_max = Some(x);
            } 
        }
    }
    
    if let (Some(x_min_c), Some(x_max_c)) = (x_min, x_max) {
        Ok((x_min_c, x_max_c))
    } else {
        Err(TextError::new("I can't compute the min and mad of an empty array"
                           .to_string()))
    }
}


#[derive(Debug, Clone)]
struct TextError { message: String }

impl TextError {
    fn new(s: String) -> TextError {
        TextError { message: s }
    }
    fn new_box(s: String) -> Box<TextError> {
        Box::new(Self::new(s))
    }
}

impl std::fmt::Display for TextError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", &self.message)
    }
}

impl std::error::Error for TextError {} 
