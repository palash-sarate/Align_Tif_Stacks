# Import required libraries
# import cv2
import easyocr
# import numpy as np
import pandas as pd
import requests
# import base64
# import gradio as gr
from gradio_client import utils as client_utils
import time

# start the easy-ocr server here and then run the below code
# docker run -it --gpus all -p 7860:7860 --platform=linux/amd64 registry.hf.space/tomofi-easyocr:latest python app.py

# Function to process the uploaded image and extract text
def ocr_from_file(image_path, roi):
    push_url = 'http://localhost:7860/api/queue/push'
    status_url = 'http://localhost:7860/api/queue/status/'
    img_b64 = client_utils.encode_url_or_file_to_base64(image_path)
    # print(img_b64[:100])
    # # print last 100 characters
    # print(img_b64[-100:])
    # Define the session hash
    session_hash = 'l1s3akiq8l9'
    # Prepare the payload
    payload = {
        "data": [img_b64, ["en"]], "action":"predict", "session_hash":session_hash
    }

    # Send the request
    r = requests.post(url=push_url, json=payload)
    # Print the response
    # print(r.json())
    r.json()
    hash = r.json().get('hash')
    # print(hash)
    # Check the status of the task

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
    
    return ocr_results

# start_time = time.time()
# ocr_from_file("C:\\Users\\Palash\\Desktop\\1_00002_med2.tif",[])
# print(f"Time taken: {time.time() - start_time} seconds")

# start_time = time.time()
# ocr_from_file("C:\\Users\\Palash\\Desktop\\1_00002_med.tif",[])
# print(f"Time taken: {time.time() - start_time} seconds")

start_time = time.time()
print(ocr_from_file("C:\\Users\\Palash\\Desktop\\1_00002.tif",[]))
print(f"Time taken: {time.time() - start_time} seconds")