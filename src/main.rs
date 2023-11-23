extern crate image;
extern crate rand;

use image::{ImageBuffer, Rgba, RgbaImage};
use rand::Rng;
use std::f64::consts::PI;
use nalgebra::Complex;


fn random_blue_to_purple_color() -> Rgba<u8> {
    let mut rng = rand::thread_rng();
    let red = rng.gen_range(0.0..0.7) * 255.0; // Up to 70% red
    let blue = rng.gen_range(0.5..1.0) * 255.0; // At least 50% blue
    Rgba([red as u8, 0, blue as u8, 255]) // Alpha is 255 for full opacity
}

fn is_in_mandelbrot_set(c: Complex<f64>, max_iterations: usize) -> bool {
    let mut z = Complex::new(0.0, 0.0);
    for _ in 0..max_iterations {
        if z.norm_sqr() > 4.0 {
            return false;
        }
        z = z * z + c;
    }
    true
}


// Constants for angles and ratios in Penrose tiling
const GOLDEN_RATIO: f64 = 1.61803398875;
const KITE_ANGLE: f64 = 2.0 * PI / 5.0; // 72 degrees
const DART_ANGLE: f64 = PI / 5.0; // 36 degrees

enum Tile {
    Kite {
        center: [f64; 2],
        size: f64,
        angle: f64,
        color: Rgba<u8>,
    },
    Dart {
        center: [f64; 2],
        size: f64,
        angle: f64,
        color: Rgba<u8>,
    },
}

fn point_at_angle_and_distance(origin: [f64; 2], angle: f64, distance: f64) -> [f64; 2] {
    [
        origin[0] + distance * angle.cos(),
        origin[1] + distance * angle.sin(),
    ]
}

fn rotate_point_about_center(point: [f64; 2], center: [f64; 2], angle: f64) -> [f64; 2] {
    // Translate point to origin
    let translated_point = [point[0] - center[0], point[1] - center[1]];

    // Rotate around origin
    let cos_theta = angle.cos();
    let sin_theta = angle.sin();
    let rotated_point = [
        cos_theta * translated_point[0] - sin_theta * translated_point[1],
        sin_theta * translated_point[0] + cos_theta * translated_point[1],
    ];

    // Translate back to original center
    [rotated_point[0] + center[0], rotated_point[1] + center[1]]
}


fn distance_from_center(point: [f64; 2], center: [f64; 2]) -> f64 {
    let dx = point[0] - center[0];
    let dy = point[1] - center[1];
    (dx * dx + dy * dy).sqrt()
}

fn interpolate_color(distance: f64, max_distance: f64) -> Rgba<u8> {
    let factor = (distance / max_distance).min(1.0);
    let red = (1.0 - factor) * 255.0; // More red for closer to the center
    let blue = (0.5 + 0.5 * factor) * 255.0; // More blue for farther from the center
    Rgba([red as u8, 0, blue as u8, 255])
}

impl Tile {
    fn to_svg_polygon(&self) -> String {
        let vertices = self.vertices();
        let points: String = vertices
            .iter()
            .map(|&[x, y]| format!("{},{}", x, y))
            .collect::<Vec<String>>()
            .join(" ");
        let color = match self {
            Tile::Kite { color, .. } | Tile::Dart { color, .. } => format!("rgb({}, {}, {})", color[0], color[1], color[2]),
        };
        format!("<polygon points=\"{}\" fill=\"{}\"/>", points, color)
    }
    // Add a method to calculate color based on distance from the center
    fn mandelbrot_color_adjustment(&self, scale: f64, offset: Complex<f64>, max_iterations: usize, max_distance: f64) -> Rgba<u8> {
        let center = self.center();
        let c = Complex::new(center[0] * scale + offset.re, center[1] * scale + offset.im);

        if is_in_mandelbrot_set(c, max_iterations) {
            // Change color if in Mandelbrot set
            Rgba([0, 255, 0, 255]) // Example: turn the tile green
        } else {
            // Use existing color calculation
            self.calculate_color(center, max_distance)
        }
    }

    fn calculate_color(&self, img_center: [f64; 2], max_distance: f64) -> Rgba<u8> {
        let distance = distance_from_center(self.center(), img_center);
        interpolate_color(distance, max_distance)
    }

    // Add a method to get the center of a tile
    fn center(&self) -> [f64; 2] {
        match *self {
            Tile::Kite { center, .. } | Tile::Dart { center, .. } => center,
        }
    }

    fn vertices(&self) -> Vec<[f64; 2]> {
        match *self {
            Tile::Kite {
                center,
                size,
                angle,
                color: _color,
            } => {
                let angle1 = angle;
                let angle2 = angle + KITE_ANGLE;
                let angle3 = angle2 + DART_ANGLE;
                let angle4 = angle3 + KITE_ANGLE;

                let point1 = point_at_angle_and_distance(center, angle1, size);
                let point2 = point_at_angle_and_distance(center, angle2, size);
                let point3 = point_at_angle_and_distance(center, angle3, size / GOLDEN_RATIO);
                let point4 = point_at_angle_and_distance(center, angle4, size / GOLDEN_RATIO);

                vec![center, point1, point2, point3, point4]
            }
            Tile::Dart {
                center,
                size,
                angle,
                color: _color,
            } => {
                let angle1 = angle;
                let angle2 = angle + DART_ANGLE;
                let angle3 = angle2 + KITE_ANGLE;
                let angle4 = angle3 + DART_ANGLE;

                let point1 = point_at_angle_and_distance(center, angle1, size);
                let point2 = point_at_angle_and_distance(center, angle2, size * GOLDEN_RATIO);
                let point3 = point_at_angle_and_distance(center, angle3, size * GOLDEN_RATIO);
                let point4 = point_at_angle_and_distance(center, angle4, size);

                vec![center, point1, point2, point3, point4]
            }
        }
    }

    fn subdivide(&self, img_center: [f64; 2], max_distance: f64) -> Vec<Tile> {
        match *self {
            Tile::Kite {
                center,
                size,
                angle,
                color: _color,
            } => {
                let new_size = size / GOLDEN_RATIO;
                let rotation_angle = PI*2.0;

                // Calculate new centers with rotation
                // Corrected function calls
                let new_center1 = rotate_point_about_center(
                    point_at_angle_and_distance(center, angle, size / 2.0),
                    center,
                    rotation_angle
                );
                let new_center2 = rotate_point_about_center(
                    point_at_angle_and_distance(center, angle + PI, size / 2.0),
                    center,
                    rotation_angle
                );
                let dart_center1 = rotate_point_about_center(
                    point_at_angle_and_distance(center, angle, size / 2.0),
                    center,
                    rotation_angle
                );
                let dart_center2 = rotate_point_about_center(
                    point_at_angle_and_distance(center, angle + PI, size / 2.0),
                    center,
                    rotation_angle
                );

                // Calculate new centers and angles
                // Adjust angles based on rotation
                let new_angle1 = angle + rotation_angle;
                let new_angle2 = angle + KITE_ANGLE + rotation_angle;
                let dart_angle1 = angle + rotation_angle;
                let dart_angle2 = angle + KITE_ANGLE + rotation_angle;

                vec![
                    Tile::Kite { center: new_center1, size: new_size, angle: new_angle1,
                        color: random_blue_to_purple_color(),
                    },
                    Tile::Kite {
                        center: new_center2,
                        size: new_size,
                        angle: new_angle2,
                        color: random_blue_to_purple_color(),
                    },
                    Tile::Kite { center: new_center2, size: new_size, angle: new_angle2,
                        color: random_blue_to_purple_color(),
                    },
                    Tile::Dart { center: dart_center1, size: new_size, angle: dart_angle1,
                        color: random_blue_to_purple_color(),
                    },
                    Tile::Dart { center: dart_center2, size: new_size, angle: dart_angle2,
                        color: random_blue_to_purple_color(),
                    },
                ]
            }
            Tile::Dart {
                center,
                size,
                angle,
                color: _color,
            } => {
                let new_size = size / GOLDEN_RATIO;

                // Calculate new centers and angles for Dart and Kite
                let new_dart_center = point_at_angle_and_distance(center, angle, size / 2.0);
                let new_dart_angle = angle;

                let new_kite_center = point_at_angle_and_distance(center, angle + PI, size / 2.0);
                let new_kite_angle = angle + DART_ANGLE;

                vec![
                    Tile::Dart {
                        center: new_dart_center,
                        size: new_size,
                        angle: new_dart_angle,
                        color: self.calculate_color(img_center, max_distance),
                    },
                    Tile::Kite {
                        center: new_kite_center,
                        size: new_size,
                        angle: new_kite_angle,
                        color: self.calculate_color(img_center, max_distance),
                    },
                ]
            }
        }
    }
}

fn subdivide_recursively(
    tile: Tile,
    depth: usize,
    img_center: [f64; 2],
    max_distance: f64,
) -> Vec<Tile> {
    if depth == 0 {
        vec![tile]
    } else {
        let mut result = Vec::new();
        for subtile in tile.subdivide(img_center, max_distance) {
            result.append(&mut subdivide_recursively(
                subtile,
                depth - 1,
                img_center,
                max_distance,
            ));
        }
        result
    }
}

fn main() {
    // Adjusted for 4K resolution
    let img_width: u32 = 3840;
    let img_height: u32 = 2160;
    let scale_factor = 2.0; // Increase this to scale the size of the tiles

    let base_center = [img_width as f64 / 2.0, img_height as f64 / 2.0];
    let offset_range = 9.0 * scale_factor; // Adjusted for scale

    let mut rng = rand::thread_rng();

    let img_center = [img_width as f64 / 2.0, img_height as f64 / 2.0];
    let max_distance = ((img_width * img_width + img_height * img_height) as f64).sqrt() / 2.0;

    let random_offset =
        |rng: &mut rand::rngs::ThreadRng| -> f64 { rng.gen_range(-offset_range..offset_range) };

    let random_angle = |rng: &mut rand::rngs::ThreadRng| -> f64 { rng.gen_range(0.0..69.0 * PI) };

    let scale_factor = 4.0; // Adjust this to scale the size of the initial tiles

    let initial_tile_size = 100.0 * scale_factor; // Increase initial tile size

    // Adjust the positions of the initial tiles
    let img_center = [img_width as f64 / 2.0, img_height as f64 / 2.0];
    let max_distance = ((img_width.pow(2) + img_height.pow(2)) as f64).sqrt() / 2.0;

    // TODO - Add Mandelbrot set parameters here


    let initial_tiles = vec![
        Tile::Kite {
            center: [
                base_center[0] + random_offset(&mut rng),
                base_center[1] + random_offset(&mut rng),
            ],
            size: initial_tile_size,
            angle: random_angle(&mut rng),
            color: interpolate_color(
                distance_from_center(
                    [
                        base_center[0] + random_offset(&mut rng),
                        base_center[1] + random_offset(&mut rng),
                    ],
                    img_center,
                ),
                max_distance,
            ),
        },
        Tile::Dart {
            center: [
                base_center[0] + initial_tile_size * GOLDEN_RATIO.cos() + random_offset(&mut rng),
                base_center[1] + initial_tile_size * GOLDEN_RATIO.sin() + random_offset(&mut rng),
            ],
            size: initial_tile_size,
            angle: KITE_ANGLE + rng.gen_range(-0.1..0.1),
            color: interpolate_color(
                distance_from_center(
                    [
                        base_center[0]
                            + initial_tile_size * GOLDEN_RATIO.cos()
                            + random_offset(&mut rng),
                        base_center[1]
                            + initial_tile_size * GOLDEN_RATIO.sin()
                            + random_offset(&mut rng),
                    ],
                    img_center,
                ),
                max_distance,
            ),
        },
        // New Kite tile
        Tile::Kite {
            center: [
                base_center[0] - 100.0 * GOLDEN_RATIO.cos() + random_offset(&mut rng),
                base_center[1] - 100.0 * GOLDEN_RATIO.sin() + random_offset(&mut rng),
            ],
            size: 100.0 * scale_factor,
            angle: -random_angle(&mut rng),
            color: interpolate_color( // TODO add Mandelbrot color adjustment here
                distance_from_center(
                    [
                        base_center[0] - 100.0 * GOLDEN_RATIO.cos() + random_offset(&mut rng),
                        base_center[1] - 100.0 * GOLDEN_RATIO.sin() + random_offset(&mut rng),
                    ],
                    img_center,
                ),
                max_distance,
            ),
        },
        // New Dart tile
        Tile::Dart {
            center: [
                base_center[0] - 50.0 + random_offset(&mut rng),
                base_center[1] + 50.0 + random_offset(&mut rng),
            ],
            size: 100.0 * scale_factor,
            angle: PI - (KITE_ANGLE + rng.gen_range(-0.1..0.1)),
            color: interpolate_color(
                distance_from_center(
                    [
                        base_center[0] - 50.0 + random_offset(&mut rng),
                        base_center[1] + 50.0 + random_offset(&mut rng),
                    ],
                    img_center,
                ),
                max_distance,
            ),
        },
        // ... Add other tiles with similar randomization if needed
    ];

    // Generate tiles recursively
    let mut tiles = Vec::new();
    for tile in initial_tiles {
        tiles.append(&mut subdivide_recursively(
            tile,
            8, // Adjust the depth if necessary
            img_center,
            max_distance,
        ));
    }
// Start building SVG content
    let mut svg_content = String::new();
    for tile in tiles {
        svg_content.push_str(&tile.to_svg_polygon());
    }

    // Create an SVG document
    let svg_data = format!("<svg width=\"{}\" height=\"{}\" xmlns=\"http://www.w3.org/2000/svg\">{}</svg>", img_width, img_height, svg_content);

    // Save to file
    std::fs::write("penrose_tiling.svg", svg_data).unwrap();
}

// Implement a function to draw polygons on the image buffer
fn draw_line(img: &mut RgbaImage, mut x0: i32, mut y0: i32, x1: i32, y1: i32, color: Rgba<u8>) {
    let dx = (x1 - x0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let dy = -((y1 - y0).abs());
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy; // error value e_xy

    loop {
        // Check bounds before drawing
        if x0 >= 0 && x0 < img.width() as i32 && y0 >= 0 && y0 < img.height() as i32 {
            img.put_pixel(x0 as u32, y0 as u32, color);
        }
        if x0 == x1 && y0 == y1 {
            break;
        }
        img.put_pixel(x0 as u32, y0 as u32, color);
        if x0 == x1 && y0 == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x0 += sx;
        } // e_xy+e_x > 0
        if e2 <= dx {
            err += dx;
            y0 += sy;
        } // e_xy+e_y < 0
    }
}

fn draw_polygon_on_image(img: &mut RgbaImage, vertices: &[[f64; 2]], color: Rgba<u8>) {
    // Draw the edges using Bresenham's algorithm
    let n = vertices.len();
    for i in 0..n {
        let (x0, y0) = (vertices[i][0] as i32, vertices[i][1] as i32);
        let (x1, y1) = (
            vertices[(i + 1) % n][0] as i32,
            vertices[(i + 1) % n][1] as i32,
        );
        draw_line(img, x0, y0, x1, y1, color);
    }

    // Polygon filling
    fill_polygon(img, vertices, color);
}

fn fill_polygon(img: &mut RgbaImage, vertices: &[[f64; 2]], color: Rgba<u8>) {
    let img_width = img.width() as f64;
    let img_height = img.height() as f64;

    // Clamp bounding box to image dimensions
    let (mut min_x, mut min_y, mut max_x, mut max_y) = (
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for &vertex in vertices {
        min_x = min_x.min(vertex[0]).max(0.0);
        min_y = min_y.min(vertex[1]).max(0.0);
        max_x = max_x.max(vertex[0]).min(img_width);
        max_y = max_y.max(vertex[1]).min(img_height);
    }

    for y in min_y as u32..=max_y as u32 {
        let mut node_x: Vec<f64> = Vec::new();

        let vertices_count = vertices.len();
        let mut j = vertices_count - 1;

        for i in 0..vertices_count {
            let (x0, y0) = (vertices[i][0], vertices[i][1]);
            let (x1, y1) = (vertices[j][0], vertices[j][1]);

            if y0 != y1 {
                let x = x0 + (y as f64 - y0) / (y1 - y0) * (x1 - x0);
                if (y0 <= y as f64 && (y as f64) < y1) || (y1 <= y as f64 && (y as f64) < y0) {
                    node_x.push(x);
                }
            }

            j = i;
        }

        // Sort the intersection points
        node_x.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Fill between pairs of intersections
        for i in (0..node_x.len()).step_by(2) {
            if i + 1 < node_x.len() {
                let start = node_x[i].max(min_x) as u32;
                let end = node_x[i + 1].min(max_x) as u32;

                for x in start..=end {
                    img.put_pixel(x, y, color);
                }
            }
        }
    }
}
