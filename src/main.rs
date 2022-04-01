use std::{env, error::Error};
use bson::Document;
use mongodb::{Client, Collection};
use mongodb::options::ClientOptions;
use serde::{Serialize, Deserialize};
use csv;
use serde_json::json;

use tokio;

#[derive(Debug, Serialize, Deserialize)]
pub struct User {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub id: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub email: Option<String>
}


fn get_users_from_csv(path: String) -> Result<Vec<User>, Box<dyn Error>> {
    let mut users: Vec<User> = vec![];
    let mut rdr = csv::Reader::from_path(path)?;
    for result in rdr.deserialize() {
        let user: User = result?;
        users.push(user);
    }
    Ok(users)
}

#[tokio::main]
async fn main() {
    dotenv::dotenv().ok();

    let database_name = env::var("DATABASE_NAME").unwrap();
    let password = env::var("PASSWORD").unwrap();
    let user = env::var("USER").unwrap();
    let file_name = env::var("FILE_NAME").unwrap();

    let host_name = format!("mongodb+srv://{}:{}@cluster0.17s4f.mongodb.net/retryWrites=true&w=majority",
                            &user,
                            &password,
    );

    let options = ClientOptions::parse(&host_name).await.unwrap();
    let client = Client::with_options(options).unwrap();
    let db = client.database(&database_name);

    let users = get_users_from_csv(file_name).unwrap();
    let mut docs: Vec<Document> = Vec::new();

    for user in users {
        let serialized_user = bson::to_bson(&user).unwrap();  // Serialize
        if let bson::Bson::Document(document) = serialized_user {
            docs.push(document)
        }
    }

    let user_collection: Collection<Document> = db.collection("users");
    let result =  user_collection.insert_many(docs, None).await.unwrap();

    println!("{}", json!(result))
}
