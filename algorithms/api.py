# pip install fastapi uvicorn[standard]
# pip install pydantic

import uvicorn
from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware

class Sequencia(BaseModel):
  seq: str

app = FastAPI()

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/analise")
async def root(response_model=Sequencia):
  return {"message":"Servidor Funcionando!"}

@app.post("/analise", response_model=Sequencia)
async def root(sequencia: Sequencia):
  return {"seq":sequencia.seq}

uvicorn.run(app, host="0.0.0.0", port=3030)