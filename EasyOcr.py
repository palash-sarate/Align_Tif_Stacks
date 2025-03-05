# Import required libraries
import cv2
# import easyocr
import numpy as np
import requests
import base64
import gradio as gr
from gradio_client import utils as client_utils
import time

# Function to process the uploaded image and extract text
def ocr_from_file(image_path, roi):
    # Read the image with OpenCV
    cropped = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    print(f'dtype: {img.dtype}, shape: {img.shape}, min: {np.min(img)}, max: {np.max(img)}')
    # Get the region of interest (ROI) from the image
    # x, y, w, h = map(int, roi)
    # cropped = cropped[y:y+h, x:x+w]
    # cropped = cv2.normalize(cropped, dst=None, alpha=0, beta=65535, norm_type=cv2.NORM_MINMAX)
    # save cropped image to disk
    # cv2.imwrite('cropped.png', cropped)
    
    # view the image to be processed
    # cv2.imshow('image', cropped)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()
    
    # Initialize the EasyOCR reader
    # reader = easyocr.Reader(['en'], gpu=True)

    # # # Perform text detection
    # results = reader.readtext(image_path)
    # # print(results)
    # # Draw bounding boxes and overlay text on the image
    # conf_threshold = 0.2
    # extracted_text = []  # To store extracted text
    # for (bbox, text, conf) in results:
    # #     # if conf > conf_threshold:
    # # #         # Append text to the list
    #     extracted_text.append(text)
    #     print(bbox, text, conf)

    #         # Get coordinates
    #         top_left = tuple(map(int, bbox[0]))
    #         bottom_right = tuple(map(int, bbox[2]))

    #         # Draw rectangle and text
    #         img = cv2.rectangle(img, top_left, bottom_right, (0, 0, 255), 2)
    #         img = cv2.putText(img, text, top_left, cv2.FONT_HERSHEY_SIMPLEX, 1, (255, 0, 0), 2)

    # Convert the image to RGB (Gradio requires RGB format)
    # img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)

    # Join the extracted text into a single string with line breaks
    # extracted_text_str = "\n".join(extracted_text)
    # print(extracted_text_str)

    # return extracted_text_str


# ocr_from_file('E:\\shake_table_data\\time_control\\1\\1_00002.tif', [0, 400, 800, 800])
# ocr_from_file('E:\\shake_table_data\\time_control\\1\\1_00002.tif', [300, 700, 100, 100])

img_path = 'E:\\shake_table_data\\time_control\\1\\1_00002.tif'
# img = cv2.imread(img_path)
# _, buffer = cv2.imencode('.png', img)
# img_b64 = base64.b64encode(buffer).decode('utf-8')
img_b64 = client_utils.encode_url_or_file_to_base64(img_path)
print(img_b64[:100])
# print last 100 characters
print(img_b64[-100:])
# Define the session hash
session_hash = 'l1s3akiq8l9'
# Prepare the payload
payload = {
    "data": [img_b64, ["en"]], "action":"predict", "session_hash":session_hash
}

# Send the request
r = requests.post(url='http://localhost:7860/api/queue/push/', json=payload)
# Print the response
print(r.json())
r.json()
hash = r.json().get('hash')
print(hash)
# Check the status of the task
status_url = 'http://localhost:7860/api/queue/status/'
while True:
    payload = {"hash": hash}
    status_response = requests.post(url = status_url, json = payload)
    status_data = status_response.json()
    if status_data.get('status') == 'COMPLETE':
        break
    elif status_data.get('status') == 'PENDING':
        print('Task is still pending...')
    else:
        print('Task failed...')
        break
    time.sleep(1)

response_data = status_data["data"]["data"]
ocr_results = response_data[1]["data"] 
print(ocr_results)
