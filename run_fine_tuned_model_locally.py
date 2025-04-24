from pathlib import Path
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM

model_dir = Path("gemma3_12b_lora_merged")   
tokenizer = AutoTokenizer.from_pretrained(model_dir, trust_remote_code=True)
model = AutoModelForCausalLM.from_pretrained(
    model_dir,
    torch_dtype=torch.float16,           # needs ≈ 24 GB VRAM; drop to int8 if <24 GB
    device_map="auto",                   # will spread across available GPUs
    load_in_8bit=True,                   # comment it if you have more vram

)

prompt = ""
inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
outputs = model.generate(**inputs, max_new_tokens=150, temperature=0.7)
print(tokenizer.decode(outputs[0], skip_special_tokens=True))