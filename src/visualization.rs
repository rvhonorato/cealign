use plotters::prelude::*;

/// Plot the alignment path, mirroring Figure 2 from Shindyalov & Bourne (1998).
///
/// Each AFP is drawn as a thick diagonal segment from (i_a, i_b) to
/// (i_a+win_size-1, i_b+win_size-1). Consecutive AFPs are joined by a thin
/// grey connector. Axes are residue indices into protein A (x) and B (y).
/// Output is written to `plot.png`.
pub fn plot(
    path: &[(usize, usize)],
    win_size: usize,
    len_a: usize,
    len_b: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("plot.png", (640, 640)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(30, 20, 20, 30);

    let mut chart = ChartBuilder::on(&root)
        .caption("Alignment path", ("sans-serif", 18))
        .x_label_area_size(30)
        .y_label_area_size(35)
        .build_cartesian_2d(0u32..len_a as u32, 0u32..len_b as u32)?;

    chart
        .configure_mesh()
        .x_desc("Protein A residue")
        .y_desc("Protein B residue")
        .x_labels(10)
        .y_labels(10)
        .x_label_formatter(&|x| format!("{}", x))
        .y_label_formatter(&|y| format!("{}", y))
        .draw()?;

    // Thin grey connectors between the end of one AFP and the start of the next
    if path.len() > 1 {
        let connectors: Vec<(u32, u32)> = path
            .windows(2)
            .flat_map(|w| {
                let (i_a, i_b) = w[0];
                let (j_a, j_b) = w[1];
                vec![
                    ((i_a + win_size - 1) as u32, (i_b + win_size - 1) as u32),
                    (j_a as u32, j_b as u32),
                ]
            })
            .collect();
        chart.draw_series(LineSeries::new(connectors, &BLACK.mix(0.35)))?;
    }

    // Thick diagonal segment for each AFP
    for &(i_a, i_b) in path {
        let segment: Vec<(u32, u32)> = (0..win_size)
            .map(|k| ((i_a + k) as u32, (i_b + k) as u32))
            .collect();
        chart.draw_series(LineSeries::new(
            segment,
            ShapeStyle::from(&BLACK).stroke_width(2),
        ))?;
    }

    root.present()?;
    Ok(())
}
