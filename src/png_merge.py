import argparse
from pathlib import Path
from PIL import Image

def merge_images(image_paths, output_path, mode="vertical"):
    if not image_paths:
        raise ValueError("At least one image path must be provided.")

    images = [Image.open(p).convert("RGBA") for p in image_paths]
    
    if mode == "vertical":
        canvas_w = max(img.width for img in images)
        canvas_h = sum(img.height for img in images)
    else:  # horizontal
        canvas_w = sum(img.width for img in images)
        canvas_h = max(img.height for img in images)

    # Create canvas with transparency
    merged_image = Image.new('RGBA', (canvas_w, canvas_h), (0, 0, 0, 0))

    offset = 0
    for img in images:
        if mode == "vertical":
            x = (canvas_w - img.width) // 2
            y = offset
            offset += img.height
        else:  # horizontal
            x = offset
            y = (canvas_h - img.height) // 2
            offset += img.width

        # Paste using image alpha channel as mask to preserve transparency
        merged_image.paste(img, (x, y), mask=img)

    # Save output
    fp = Path(output_path)
    fp.parent.mkdir(parents=True, exist_ok=True)
    merged_image.save(output_path)
    print(f"Successfully saved {mode} merged image ({len(images)} files) to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multiple PNG images vertically or horizontally.")
    parser.add_argument(
        '--images',
        required=True,
        type=lambda s: [p.strip() for p in s.split(',') if p.strip()],
        help="Comma-separated list of image paths (e.g., 'a.png,b.png,c.png')"
    )
    parser.add_argument('--output', required=True, help="Path for the output merged image")
    parser.add_argument(
        '--mode', 
        choices=['vertical', 'horizontal'], 
        default='vertical', 
        help="Merge orientation: 'vertical' (default) or 'horizontal'"
    )
    
    args = parser.parse_args()
    merge_images(args.images, args.output, args.mode)
